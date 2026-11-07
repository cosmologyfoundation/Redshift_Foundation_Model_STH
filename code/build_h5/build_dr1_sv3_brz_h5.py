"""Build a filtered DESI DR1 SV3 H5 dataset with coadded BRZ spectra.

The catalog selection is exactly:

    SURVEY == "sv3"
    ZCAT_PRIMARY == True
    OBJTYPE == "TGT"
    ZWARN == 0

Selected rows are grouped by SURVEY, PROGRAM, and HEALPIX so that each coadd
FITS file is opened once. B/R/Z cameras are combined with
desispec.coadd_cameras, and the output H5 stores BRZ_FLUX, BRZ_IVAR, and
BRZ_MASK without duplicating the shared wavelength grid.
"""

import argparse
import hashlib
import os
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

import h5py
import numpy as np
from astropy.io import fits
from astropy.table import Table
from desispec.coaddition import coadd_cameras
from desispec.spectra import Spectra

try:
    import desispec

    DESISPEC_VERSION = getattr(desispec, "__version__", "unknown")
except Exception:
    DESISPEC_VERSION = "unknown"


DEFAULT_ZALL = Path(
    "/global/cfs/cdirs/desi/public/dr1/spectro/redux/iron/"
    "zcatalog/v1/zall-pix-iron.fits"
)
DEFAULT_COADD_ROOT = Path(
    "/global/cfs/cdirs/desi/public/dr1/spectro/redux/iron/healpix"
)
DEFAULT_OUTPUT = Path(
    "/pscratch/sd/s/shunth/dr1_sv3_primary_tgt_zwarn0_brz.h5"
)

FILTER_DESCRIPTION = (
    'SURVEY == "sv3" and ZCAT_PRIMARY == True and '
    'OBJTYPE == "TGT" and ZWARN == 0'
)

ZALL_FIELDS = (
    "TARGETID",
    "SURVEY",
    "PROGRAM",
    "HEALPIX",
    "Z",
    "ZERR",
    "SPECTYPE",
    "SUBTYPE",
    "ZWARN",
    "CHI2",
    "DELTACHI2",
    "NPIXELS",
    "COADD_FIBERSTATUS",
    "OBJTYPE",
    "SV3_DESI_TARGET",
    "SV3_BGS_TARGET",
    "SV3_MWS_TARGET",
    "SV3_SCND_TARGET",
    "ZCAT_PRIMARY",
    "ZCAT_NSPEC",
    "SV_PRIMARY",
    "SV_NSPEC",
)

STRING_FIELDS = {"SURVEY", "PROGRAM", "SPECTYPE", "SUBTYPE", "OBJTYPE"}
BOOL_FIELDS = {"ZCAT_PRIMARY", "SV_PRIMARY"}

# Preserve catalog precision for redshift labels and quality statistics.
H5_DTYPES = {
    "TARGETID": np.int64,
    "HEALPIX": np.int32,
    "Z": np.float64,
    "ZERR": np.float64,
    "ZWARN": np.int64,
    "CHI2": np.float64,
    "DELTACHI2": np.float64,
    "NPIXELS": np.int64,
    "COADD_FIBERSTATUS": np.int64,
    "SV3_DESI_TARGET": np.int64,
    "SV3_BGS_TARGET": np.int64,
    "SV3_MWS_TARGET": np.int64,
    "SV3_SCND_TARGET": np.int64,
    "ZCAT_PRIMARY": np.bool_,
    "ZCAT_NSPEC": np.int32,
    "SV_PRIMARY": np.bool_,
    "SV_NSPEC": np.int32,
}

SPECTRUM_FIELDS = {
    "BRZ_FLUX": np.float32,
    "BRZ_IVAR": np.float32,
    "BRZ_MASK": np.uint32,
}


def clean_strings(values: np.ndarray) -> np.ndarray:
    """Decode and trim a FITS string column into a NumPy unicode array."""
    return np.char.strip(np.asarray(values).astype(str))


def select_catalog_rows(
    zcat: np.ndarray,
    chunk_size: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return selected catalog indices and their grouping columns."""
    index_parts: List[np.ndarray] = []
    survey_parts: List[np.ndarray] = []
    program_parts: List[np.ndarray] = []
    healpix_parts: List[np.ndarray] = []
    total = len(zcat)

    for start in range(0, total, chunk_size):
        stop = min(start + chunk_size, total)
        survey = clean_strings(zcat["SURVEY"][start:stop])
        objtype = clean_strings(zcat["OBJTYPE"][start:stop])
        good = (
            (survey == "sv3")
            & np.asarray(zcat["ZCAT_PRIMARY"][start:stop], dtype=bool)
            & (objtype == "TGT")
            & (np.asarray(zcat["ZWARN"][start:stop]) == 0)
        )

        local = np.flatnonzero(good)
        if local.size:
            absolute = local.astype(np.int64, copy=False) + start
            index_parts.append(absolute)
            # Compact fixed-width bytes keep the 20M-row grouping arrays small.
            survey_parts.append(np.asarray(survey[local], dtype="S7"))
            program_parts.append(
                np.asarray(clean_strings(zcat["PROGRAM"][absolute]), dtype="S6")
            )
            healpix_parts.append(
                np.asarray(zcat["HEALPIX"][absolute], dtype=np.int32)
            )

        print(
            f"catalog scan: {stop:,}/{total:,} | "
            f"selected={sum(part.size for part in index_parts):,}",
            flush=True,
        )

    if not index_parts:
        raise RuntimeError("The requested catalog selection returned zero rows")

    indices = np.concatenate(index_parts)
    surveys = np.concatenate(survey_parts)
    programs = np.concatenate(program_parts)
    healpix = np.concatenate(healpix_parts)

    # Stable catalog-index tie-breaker makes resume deterministic.
    order = np.lexsort((indices, healpix, programs, surveys))
    return indices[order], surveys[order], programs[order], healpix[order]


def iter_groups(
    indices: np.ndarray,
    surveys: np.ndarray,
    programs: np.ndarray,
    healpix: np.ndarray,
) -> Iterator[Tuple[str, str, int, np.ndarray]]:
    """Yield contiguous SURVEY/PROGRAM/HEALPIX groups."""
    changes = np.flatnonzero(
        (surveys[1:] != surveys[:-1])
        | (programs[1:] != programs[:-1])
        | (healpix[1:] != healpix[:-1])
    ) + 1
    starts = np.concatenate(([0], changes))
    stops = np.concatenate((changes, [indices.size]))

    for start, stop in zip(starts, stops):
        survey = surveys[start].decode("ascii").strip()
        program = programs[start].decode("ascii").strip()
        yield survey, program, int(healpix[start]), indices[start:stop]


def coadd_path(root: Path, survey: str, program: str, healpix: int) -> Path:
    hpixgroup = healpix // 100
    return (
        root
        / survey
        / program
        / str(hpixgroup)
        / str(healpix)
        / f"coadd-{survey}-{program}-{healpix}.fits"
    )


def compression_options(name: str, gzip_level: int) -> Dict[str, object]:
    if name == "none":
        return {}
    if name == "lzf":
        return {"compression": "lzf", "shuffle": True}
    return {
        "compression": "gzip",
        "compression_opts": gzip_level,
        "shuffle": True,
    }


def create_output(
    path: Path,
    source_zall: Path,
    selected_rows: int,
    metadata_chunk_rows: int,
) -> h5py.File:
    path.parent.mkdir(parents=True, exist_ok=True)
    h5 = h5py.File(path, "w")
    string_dtype = h5py.string_dtype(encoding="utf-8")

    for field in ZALL_FIELDS:
        dtype = string_dtype if field in STRING_FIELDS else H5_DTYPES[field]
        h5.create_dataset(
            field,
            shape=(0,),
            maxshape=(None,),
            chunks=(metadata_chunk_rows,),
            dtype=dtype,
        )

    h5.attrs["SOURCE_ZCATALOG"] = str(source_zall)
    h5.attrs["CATALOG_FILTER"] = FILTER_DESCRIPTION
    h5.attrs["DESISPEC_VERSION"] = DESISPEC_VERSION
    h5.attrs["SELECTED_CATALOG_ROWS"] = selected_rows
    h5.attrs["ROWS_COMMITTED"] = 0
    h5.attrs["NEXT_GROUP_INDEX"] = 0
    h5.attrs["BUILD_COMPLETE"] = False
    h5.flush()
    return h5


def ensure_spectrum_datasets(
    h5: h5py.File,
    npix: int,
    chunk_rows: int,
    compression: str,
    gzip_level: int,
) -> None:
    options = compression_options(compression, gzip_level)
    for field, dtype in SPECTRUM_FIELDS.items():
        if field in h5:
            if h5[field].shape[1] != npix:
                raise RuntimeError(
                    f"Existing {field} has {h5[field].shape[1]} pixels, expected {npix}"
                )
            continue
        h5.create_dataset(
            field,
            shape=(0, npix),
            maxshape=(None, npix),
            chunks=(chunk_rows, npix),
            dtype=dtype,
            **options,
        )


def truncate_to_committed_rows(h5: h5py.File) -> int:
    committed = int(h5.attrs.get("ROWS_COMMITTED", 0))
    for dataset in h5.values():
        if dataset.shape[0] != committed:
            dataset.resize(committed, axis=0)
    h5.flush()
    return committed


def append_batch(
    h5: h5py.File,
    zcat: np.ndarray,
    catalog_indices: np.ndarray,
    combined: Spectra,
) -> None:
    batch_size = catalog_indices.size
    old_size = h5["TARGETID"].shape[0]
    new_size = old_size + batch_size

    for field in ZALL_FIELDS:
        dataset = h5[field]
        dataset.resize(new_size, axis=0)
        values = zcat[field][catalog_indices]
        if field in STRING_FIELDS:
            values = np.asarray(clean_strings(values), dtype=object)
        elif field in BOOL_FIELDS:
            values = np.asarray(values, dtype=bool)
        else:
            values = np.asarray(values, dtype=H5_DTYPES[field])
        dataset[old_size:new_size] = values

    spectra = {
        "BRZ_FLUX": np.asarray(combined.flux["brz"], dtype=np.float32),
        "BRZ_IVAR": np.asarray(combined.ivar["brz"], dtype=np.float32),
        "BRZ_MASK": np.asarray(combined.mask["brz"], dtype=np.uint32),
    }
    for field, values in spectra.items():
        dataset = h5[field]
        dataset.resize(new_size, axis=0)
        dataset[old_size:new_size] = values


def match_fibermap_rows(
    fibermap_targetids: np.ndarray,
    requested_targetids: np.ndarray,
) -> np.ndarray:
    """Map TARGETIDs to coadd rows without silently accepting duplicates."""
    fibermap_targetids = np.asarray(fibermap_targetids, dtype=np.int64)
    requested_targetids = np.asarray(requested_targetids, dtype=np.int64)
    order = np.argsort(fibermap_targetids)
    sorted_ids = fibermap_targetids[order]

    if sorted_ids.size > 1 and np.any(sorted_ids[1:] == sorted_ids[:-1]):
        raise RuntimeError("FIBERMAP contains duplicate TARGETID values")

    positions = np.searchsorted(sorted_ids, requested_targetids)
    inside = positions < sorted_ids.size
    matches = np.zeros(requested_targetids.size, dtype=bool)
    matches[inside] = sorted_ids[positions[inside]] == requested_targetids[inside]
    if not np.all(matches):
        missing = requested_targetids[~matches][:10]
        raise RuntimeError(
            f"{np.count_nonzero(~matches)} TARGETIDs are absent from FIBERMAP; "
            f"first missing values: {missing.tolist()}"
        )
    return order[positions]


def make_spectra(coadd: fits.HDUList, rows: np.ndarray) -> Spectra:
    bands = ("b", "r", "z")
    # desispec 0.68.x expects Astropy Table.colnames in coadd_cameras().
    fibermap = Table(coadd["FIBERMAP"].data[rows], copy=True)
    return Spectra(
        bands=list(bands),
        wave={
            band: np.asarray(coadd[f"{band.upper()}_WAVELENGTH"].data)
            for band in bands
        },
        flux={
            band: np.asarray(
                coadd[f"{band.upper()}_FLUX"].data[rows], dtype=np.float32
            )
            for band in bands
        },
        ivar={
            band: np.asarray(
                coadd[f"{band.upper()}_IVAR"].data[rows], dtype=np.float32
            )
            for band in bands
        },
        mask={
            band: np.asarray(
                coadd[f"{band.upper()}_MASK"].data[rows], dtype=np.uint32
            )
            for band in bands
        },
        fibermap=fibermap,
        single=True,
    )


def validate_wave(h5: h5py.File, wave: np.ndarray) -> None:
    wave = np.asarray(wave, dtype=np.float64)
    if wave.ndim != 1 or wave.size == 0 or not np.all(np.diff(wave) > 0):
        raise RuntimeError("coadd_cameras returned an invalid BRZ wavelength grid")

    digest = hashlib.sha256(wave.tobytes()).hexdigest()
    if "BRZ_WAVE_SHA256" in h5.attrs:
        if str(h5.attrs["BRZ_WAVE_SHA256"]) != digest:
            raise RuntimeError("BRZ wavelength grid changed between coadd files")
        return

    differences = np.diff(wave)
    h5.attrs["BRZ_NPIX"] = wave.size
    h5.attrs["BRZ_WAVE_MIN"] = float(wave[0])
    h5.attrs["BRZ_WAVE_MAX"] = float(wave[-1])
    h5.attrs["BRZ_WAVE_STEP"] = float(np.median(differences))
    h5.attrs["BRZ_WAVE_SHA256"] = digest


def validate_existing_output(h5: h5py.File, source_zall: Path) -> None:
    if str(h5.attrs.get("SOURCE_ZCATALOG", "")) != str(source_zall):
        raise RuntimeError("Existing output was built from a different zcatalog")
    if str(h5.attrs.get("CATALOG_FILTER", "")) != FILTER_DESCRIPTION:
        raise RuntimeError("Existing output uses a different catalog filter")
    missing = [field for field in ZALL_FIELDS if field not in h5]
    if missing:
        raise RuntimeError(f"Existing output is missing metadata fields: {missing}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zall", type=Path, default=DEFAULT_ZALL)
    parser.add_argument("--coadd-root", type=Path, default=DEFAULT_COADD_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--catalog-chunk-size", type=int, default=1_000_000)
    parser.add_argument("--coadd-batch-size", type=int, default=128)
    parser.add_argument("--h5-chunk-rows", type=int, default=64)
    parser.add_argument("--metadata-chunk-rows", type=int, default=8192)
    parser.add_argument(
        "--compression", choices=("lzf", "gzip", "none"), default="lzf"
    )
    parser.add_argument("--gzip-level", type=int, default=4)
    parser.add_argument("--max-groups", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.overwrite and args.resume:
        raise ValueError("Choose either --overwrite or --resume, not both")
    if args.catalog_chunk_size <= 0 or args.coadd_batch_size <= 0:
        raise ValueError("Batch and chunk sizes must be positive")
    if not args.zall.is_file():
        raise FileNotFoundError(f"zcatalog not found: {args.zall}")
    if not args.coadd_root.is_dir():
        raise FileNotFoundError(f"coadd root not found: {args.coadd_root}")

    with fits.open(args.zall, memmap=True, lazy_load_hdus=True) as zall_hdul:
        zcat_hdu = zall_hdul["ZCATALOG"]
        available = set(zcat_hdu.columns.names)
        missing_fields = [field for field in ZALL_FIELDS if field not in available]
        if missing_fields:
            raise RuntimeError(f"ZCATALOG is missing required fields: {missing_fields}")
        zcat = zcat_hdu.data

        selected, surveys, programs, healpix = select_catalog_rows(
            zcat, args.catalog_chunk_size
        )
        groups = list(iter_groups(selected, surveys, programs, healpix))
        unique_targetids = np.unique(
            np.asarray(zcat["TARGETID"][selected], dtype=np.int64)
        ).size
        print(f"selected rows: {selected.size:,}")
        print(f"unique TARGETIDs: {unique_targetids:,}")
        print(f"coadd groups: {len(groups):,}")
        if unique_targetids != selected.size:
            raise RuntimeError(
                "ZCAT_PRIMARY selection still contains duplicate TARGETIDs"
            )
        if args.dry_run:
            return

        if args.output.exists() and args.overwrite:
            os.remove(args.output)
        if args.output.exists() and not args.resume:
            raise FileExistsError(
                f"Output exists: {args.output}; use --resume or --overwrite"
            )

        if args.output.exists():
            h5 = h5py.File(args.output, "r+")
            validate_existing_output(h5, args.zall)
            truncate_to_committed_rows(h5)
            start_group = int(h5.attrs.get("NEXT_GROUP_INDEX", 0))
            print(
                f"resuming at group {start_group:,}, "
                f"committed rows={int(h5.attrs['ROWS_COMMITTED']):,}"
            )
        else:
            h5 = create_output(
                args.output,
                args.zall,
                selected.size,
                args.metadata_chunk_rows,
            )
            start_group = 0

        stop_group = len(groups)
        if args.max_groups is not None:
            stop_group = min(stop_group, start_group + args.max_groups)

        try:
            for group_index in range(start_group, stop_group):
                survey, program, hpix, catalog_indices = groups[group_index]
                path = coadd_path(args.coadd_root, survey, program, hpix)
                if not path.is_file():
                    raise FileNotFoundError(f"coadd file not found: {path}")

                # MASK image HDUs use unsigned-integer FITS scaling (BZERO), which
                # Astropy cannot expose through a memory map.
                with fits.open(path, memmap=False, lazy_load_hdus=True) as coadd:
                    fibermap_ids = coadd["FIBERMAP"].data["TARGETID"]
                    requested_ids = np.asarray(
                        zcat["TARGETID"][catalog_indices], dtype=np.int64
                    )
                    coadd_rows = match_fibermap_rows(fibermap_ids, requested_ids)

                    for offset in range(0, catalog_indices.size, args.coadd_batch_size):
                        stop = min(offset + args.coadd_batch_size, catalog_indices.size)
                        catalog_batch = catalog_indices[offset:stop]
                        coadd_batch = coadd_rows[offset:stop]
                        combined = coadd_cameras(make_spectra(coadd, coadd_batch))
                        expected_ids = np.asarray(
                            zcat["TARGETID"][catalog_batch], dtype=np.int64
                        )
                        combined_ids = np.asarray(
                            combined.fibermap["TARGETID"], dtype=np.int64
                        )
                        if not np.array_equal(combined_ids, expected_ids):
                            raise RuntimeError(
                                "coadd_cameras changed TARGETID order; refusing to "
                                "write misaligned spectra and catalog metadata"
                            )
                        wave = np.asarray(combined.wave["brz"])
                        validate_wave(h5, wave)
                        ensure_spectrum_datasets(
                            h5,
                            wave.size,
                            args.h5_chunk_rows,
                            args.compression,
                            args.gzip_level,
                        )
                        append_batch(h5, zcat, catalog_batch, combined)

                committed = h5["TARGETID"].shape[0]
                h5.attrs["ROWS_COMMITTED"] = committed
                h5.attrs["NEXT_GROUP_INDEX"] = group_index + 1
                h5.flush()

                if (group_index + 1) % 10 == 0 or group_index + 1 == stop_group:
                    print(
                        f"groups: {group_index + 1:,}/{len(groups):,} | "
                        f"rows: {committed:,}/{selected.size:,}",
                        flush=True,
                    )

            complete = stop_group == len(groups)
            if complete and h5["TARGETID"].shape[0] != selected.size:
                raise RuntimeError(
                    f"Completed all groups but wrote {h5['TARGETID'].shape[0]:,} "
                    f"of {selected.size:,} selected rows"
                )
            h5.attrs["BUILD_COMPLETE"] = complete
            h5.attrs["ROWS_WRITTEN"] = h5["TARGETID"].shape[0]
            h5.flush()
            print(f"output: {args.output}")
            print(f"rows written: {h5['TARGETID'].shape[0]:,}")
            print(f"build complete: {complete}")
        finally:
            h5.close()


if __name__ == "__main__":
    main()
