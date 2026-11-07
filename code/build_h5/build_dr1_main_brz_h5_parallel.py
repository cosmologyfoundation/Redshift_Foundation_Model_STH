#!/usr/bin/env python3
"""Parallel MAIN counterpart of build_dr1_sv3_brz_h5.py.

Selection: SURVEY=main, ZCAT_PRIMARY, OBJTYPE=TGT, ZWARN=0.
This is a labeled catalog export, not a Redrock-blind cohort selection.
Uses the unchanged original builder's official DESI camera coaddition and
HDF5 layout. MAIN targeting/primary metadata replace the SV3/SV columns.
CPU processes scan disjoint catalog ranges, build/ compress separate shards,
and validate disjoint output ranges. Only the final HDF5 byte-copy writer
is serialized (ordinary HDF5 must not have competing writers).
All but the last shard have a multiple of 64 rows, so compressed chunks can
be copied without serial recompression. No VDS/external shard dependency.
Resumable at catalog chunks, shard coadd groups, and merged shard boundaries.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import time

import fitsio
import h5py
import numpy as np
from astropy.io import fits
import build_dr1_sv3_brz_h5 as base


RENAMES = {"SV3_" + x: x for x in ("DESI_TARGET", "BGS_TARGET", "MWS_TARGET", "SCND_TARGET")}
RENAMES.update(SV_PRIMARY="MAIN_PRIMARY", SV_NSPEC="MAIN_NSPEC")
base.ZALL_FIELDS = tuple(RENAMES.get(k, k) for k in base.ZALL_FIELDS)
base.H5_DTYPES = {RENAMES.get(k, k): v for k, v in base.H5_DTYPES.items()}
base.BOOL_FIELDS = {RENAMES.get(k, k) for k in base.BOOL_FIELDS}
base.FILTER_DESCRIPTION = 'SURVEY == "main" and ZCAT_PRIMARY == True and OBJTYPE == "TGT" and ZWARN == 0'
CHUNK = 64
CTX = mp.get_context("spawn")


def log(message):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), message, flush=True)


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def selection_ok(a):
    return ((base.clean_strings(a["SURVEY"]) == "main")
            & a["ZCAT_PRIMARY"].astype(bool)
            & (base.clean_strings(a["OBJTYPE"]) == "TGT")
            & (a["ZWARN"] == 0))


def scan_chunk(task):
    source, start, stop, directory = task
    path = Path(directory) / f"scan_{start:09d}_{stop:09d}.npy"
    if path.exists():
        return str(path)
    with fitsio.FITS(source) as f:
        a = f["ZCATALOG"].read(columns=list(base.ZALL_FIELDS), rows=np.arange(start, stop))
    a = a[selection_ok(a)]
    temporary = path.with_suffix(".tmp.npy")
    np.save(temporary, a, allow_pickle=False)
    temporary.replace(path)
    log(f"catalog {start:,}:{stop:,} selected {len(a):,}")
    return str(path)


def group_slices(a):
    changes = np.flatnonzero((a["SURVEY"][1:] != a["SURVEY"][:-1])
                            | (a["PROGRAM"][1:] != a["PROGRAM"][:-1])
                            | (a["HEALPIX"][1:] != a["HEALPIX"][:-1])) + 1
    edges = np.r_[0, changes, len(a)]
    return list(zip(edges[:-1], edges[1:]))


def prepare_selection(args):
    directory = args.work_dir / "catalog"
    directory.mkdir(parents=True, exist_ok=True)
    config = dict(source=str(args.zall), source_size=args.zall.stat().st_size,
                  source_mtime_ns=args.zall.stat().st_mtime_ns,
                  fields=list(base.ZALL_FIELDS), filter=base.FILTER_DESCRIPTION)
    config_path = directory / "config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise RuntimeError("Catalog cache configuration mismatch; use a new work directory")
    atomic_json(config_path, config)
    path = directory / "selection.npy"
    if not path.exists():
        with fitsio.FITS(str(args.zall)) as f:
            missing = set(base.ZALL_FIELDS) - set(f["ZCATALOG"].get_colnames())
            if missing:
                raise RuntimeError(f"Missing catalog columns: {missing}")
            n = f["ZCATALOG"].get_nrows()
        tasks = [(str(args.zall), start, min(start + args.catalog_chunk_rows, n), str(directory))
                 for start in range(0, n, args.catalog_chunk_rows)]
        with ProcessPoolExecutor(max_workers=min(8, args.workers), mp_context=CTX) as pool:
            paths = list(pool.map(scan_chunk, tasks))
        a = np.concatenate([np.load(p, allow_pickle=False) for p in paths])
        if not len(a) or not np.all(selection_ok(a)):
            raise RuntimeError("Empty or invalid selection")
        if len(np.unique(a["TARGETID"])) != len(a):
            raise RuntimeError("ZCAT_PRIMARY selection contains duplicate TARGETID")
        order = np.lexsort((a["TARGETID"], a["HEALPIX"], a["PROGRAM"], a["SURVEY"]))
        a = a[order]
        np.save(directory / "selection.tmp.npy", a, allow_pickle=False)
        (directory / "selection.tmp.npy").replace(path)
        atomic_json(directory / "summary.json", dict(rows=len(a), unique_targetids=len(a),
                    groups=len(group_slices(a)), filter=base.FILTER_DESCRIPTION,
                    programs={str(k): int(v) for k, v in zip(*np.unique(a["PROGRAM"], return_counts=True))}))
    a = np.load(path, mmap_mode="r", allow_pickle=False)
    if args.smoke_groups:
        groups = group_slices(a)
        chosen = np.unique(np.linspace(0, len(groups) - 1, min(args.smoke_groups, len(groups)), dtype=int))
        subset = np.concatenate([a[groups[i][0]:groups[i][1]] for i in chosen])
        path = directory / f"smoke_{args.smoke_groups}.npy"
        if not path.exists():
            np.save(path, subset, allow_pickle=False)
        else:
            assert np.array_equal(np.load(path)["TARGETID"], subset["TARGETID"])
        a = np.load(path, mmap_mode="r", allow_pickle=False)
    log(f"selection rows={len(a):,}, coadd groups={len(group_slices(a)):,}")
    return path, a


def build_shard(task):
    manifest, start, stop, shard_name, source, root = task
    a = np.load(manifest, mmap_mode="r", allow_pickle=False)[start:stop]
    path = Path(shard_name)
    groups = group_slices(a)
    if path.exists():
        h = h5py.File(path, "r+")
        base.validate_existing_output(h, Path(source))
        if int(h.attrs["SHARD_START"]) != start or int(h.attrs["SHARD_STOP"]) != stop:
            raise RuntimeError("Shard range mismatch")
        if h.attrs["DESISPEC_VERSION"] != base.DESISPEC_VERSION:
            raise RuntimeError("Cannot resume with a different DESI version")
        if h.attrs["BUILD_COMPLETE"]:
            h.close()
            return shard_name
        base.truncate_to_committed_rows(h)
    else:
        h = base.create_output(path, Path(source), len(a), 8192)
        h.attrs.update(SHARD_START=start, SHARD_STOP=stop)
    try:
        first_group = int(h.attrs["NEXT_GROUP_INDEX"])
        for gi in range(first_group, len(groups)):
            lo, hi = groups[gi]
            survey = str(base.clean_strings(a["SURVEY"][lo:lo + 1])[0])
            program = str(base.clean_strings(a["PROGRAM"][lo:lo + 1])[0])
            p = base.coadd_path(Path(root), survey, program, int(a["HEALPIX"][lo]))
            with fits.open(p, memmap=False, lazy_load_hdus=True) as coadd:
                rows = base.match_fibermap_rows(coadd["FIBERMAP"].data["TARGETID"], a["TARGETID"][lo:hi])
                for offset in range(0, hi - lo, 128):
                    idx = np.arange(lo + offset, min(lo + offset + 128, hi))
                    combined = base.coadd_cameras(base.make_spectra(coadd, rows[offset:offset + len(idx)]))
                    if not np.array_equal(combined.fibermap["TARGETID"], a["TARGETID"][idx]):
                        raise RuntimeError("Official coaddition changed TARGETID order")
                    base.validate_wave(h, combined.wave["brz"])
                    base.ensure_spectrum_datasets(h, len(combined.wave["brz"]), CHUNK, "lzf", 4)
                    base.append_batch(h, a, idx, combined)
            # Flush data BEFORE advancing the durable resume marker.
            h.flush()
            h.attrs["ROWS_COMMITTED"] = hi
            h.attrs["NEXT_GROUP_INDEX"] = gi + 1
            h.flush()
            if (gi + 1) % 20 == 0 or gi + 1 == len(groups):
                log(f"{path.name}: groups {gi + 1}/{len(groups)}, rows {hi:,}/{len(a):,}")
        if len(h["TARGETID"]) != len(a):
            raise RuntimeError("Incomplete shard")
        h.attrs["BUILD_COMPLETE"] = True
        h.attrs["ROWS_WRITTEN"] = len(a)
        h.flush()
    finally:
        h.close()
    return shard_name


def merge_shards(args, tasks, total, groups):
    partial = args.output.with_suffix(".partial.h5")
    if args.output.exists():
        with h5py.File(args.output, "r") as h:
            if not h.attrs.get("BUILD_COMPLETE", False):
                raise RuntimeError("Existing final output is not complete")
            if h.attrs.get("BUILD_SIGNATURE", "") != args.signature:
                raise RuntimeError("Existing final output is from another configuration")
        log("Final output already complete; validating again")
        return args.output
    if partial.exists():
        h = h5py.File(partial, "r+")
        base.validate_existing_output(h, args.zall)
        if h.attrs.get("BUILD_SIGNATURE", "") != args.signature:
            raise RuntimeError("Partial merge configuration mismatch")
        base.truncate_to_committed_rows(h)
    else:
        h = base.create_output(partial, args.zall, total, 8192)
        h.attrs.update(MERGED_SHARDS=0, BUILD_SIGNATURE=args.signature,
                       CPU_WORKERS=args.workers, SHARD_COUNT=len(tasks),
                       BUILD_IMPLEMENTATION="parallel shards; raw compressed chunk merge",
                       ORIGINAL_SV3_SCRIPT_SHA256=hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),
                       BUILD_SCRIPT_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       SMOKE_TEST=bool(args.smoke_groups))
    try:
        for shard_i in range(int(h.attrs["MERGED_SHARDS"]), len(tasks)):
            _, start, stop, shard_path, _, _ = tasks[shard_i]
            with h5py.File(shard_path, "r") as s:
                if not s.attrs["BUILD_COMPLETE"] or len(s["TARGETID"]) != stop - start:
                    raise RuntimeError("Incomplete source shard")
                if start % CHUNK:
                    raise RuntimeError("Shard is not HDF5-chunk aligned")
                if "BRZ_WAVE_SHA256" in h.attrs and h.attrs["BRZ_WAVE_SHA256"] != s.attrs["BRZ_WAVE_SHA256"]:
                    raise RuntimeError("Shard wavelength grid mismatch")
                for key, value in s.attrs.items():
                    if key.startswith("BRZ_"):
                        h.attrs[key] = value
                base.ensure_spectrum_datasets(h, int(s.attrs["BRZ_NPIX"]), CHUNK, "lzf", 4)
                for field in base.ZALL_FIELDS:
                    h[field].resize(stop, axis=0)
                    h[field][start:stop] = s[field][:]
                for field in base.SPECTRUM_FIELDS:
                    dest, src = h[field], s[field]
                    dest.resize(stop, axis=0)
                    if dest.chunks != src.chunks:
                        raise RuntimeError("Chunk layout mismatch")
                    for r in range(0, stop - start, CHUNK):
                        mask, raw = src.id.read_direct_chunk((r, 0))
                        dest.id.write_direct_chunk((start + r, 0), raw, filter_mask=mask)
                h.flush()
                h.attrs["ROWS_COMMITTED"] = stop
                h.attrs["MERGED_SHARDS"] = shard_i + 1
                h.flush()
            log(f"merged shard {shard_i + 1}/{len(tasks)}, rows={stop:,}/{total:,}")
        h.attrs["NEXT_GROUP_INDEX"] = groups
        h.attrs["ROWS_WRITTEN"] = total
        h.flush()
    finally:
        h.close()
    return partial


def validate_shard(task):
    build_task, output = task
    manifest, start, stop, shard_path, _, root = build_task
    a = np.load(manifest, mmap_mode="r", allow_pickle=False)[start:stop]
    assert np.all(selection_ok(a))
    with h5py.File(output, "r") as h, h5py.File(shard_path, "r") as s:
        for field in base.ZALL_FIELDS:
            actual = h[field][start:stop]
            expected = a[field]
            if field in base.STRING_FIELDS:
                equal = np.array_equal(base.clean_strings(actual), base.clean_strings(expected))
            elif expected.dtype.kind == "f":
                equal = np.array_equal(actual, expected, equal_nan=True)
            else:
                equal = np.array_equal(actual, expected)
            if not equal:
                raise RuntimeError(f"Metadata mismatch: {field} {start}:{stop}")
        # Compare every compressed spectral chunk, bit-for-bit, after merge.
        digest = hashlib.sha256()
        for field in base.SPECTRUM_FIELDS:
            for r in range(0, stop - start, CHUNK):
                sm, sb = s[field].id.read_direct_chunk((r, 0))
                hm, hb = h[field].id.read_direct_chunk((start + r, 0))
                if sm != hm or sb != hb:
                    raise RuntimeError(f"Spectral chunk mismatch: {field} row {start+r}")
                digest.update(hb)
        # Independently re-open source FITS and recompute endpoint spectra.
        for r in sorted({0, len(a) - 1}):
            survey = str(base.clean_strings(a["SURVEY"][r:r + 1])[0])
            program = str(base.clean_strings(a["PROGRAM"][r:r + 1])[0])
            p = base.coadd_path(Path(root), survey, program, int(a["HEALPIX"][r]))
            with fits.open(p, memmap=False, lazy_load_hdus=True) as f:
                rows = base.match_fibermap_rows(f["FIBERMAP"].data["TARGETID"], a["TARGETID"][r:r + 1])
                c = base.coadd_cameras(base.make_spectra(f, rows))
                for field, values in (("BRZ_FLUX", c.flux["brz"]), ("BRZ_IVAR", c.ivar["brz"]), ("BRZ_MASK", c.mask["brz"])):
                    if not np.array_equal(h[field][start + r], values[0], equal_nan=True):
                        raise RuntimeError(f"Source recomputation mismatch {field} row {start+r}")
        return dict(start=start, stop=stop, rows=stop-start, compressed_sha256=digest.hexdigest(),
                    metadata_exact=True, all_spectral_chunks_exact=True, source_recomputed_rows=min(2,len(a)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--zall", type=Path, default=base.DEFAULT_ZALL)
    p.add_argument("--coadd-root", type=Path, default=base.DEFAULT_COADD_ROOT)
    p.add_argument("--output", type=Path, default=Path("/pscratch/sd/s/shunth/dr1_main_primary_tgt_zwarn0_brz.h5"))
    p.add_argument("--work-dir", type=Path, default=Path("/pscratch/sd/s/shunth/codex_learning_map/dr1_main_brz_build"))
    p.add_argument("--workers", type=int, default=64)
    p.add_argument("--catalog-chunk-rows", type=int, default=2_000_000)
    p.add_argument("--smoke-groups", type=int, default=0)
    p.add_argument("--prepare-only", action="store_true")
    args = p.parse_args()
    if args.workers < 2 or args.catalog_chunk_rows <= 0 or args.smoke_groups < 0:
        p.error("At least 2 worker processes and positive chunk size are required")
    if args.workers > len(os.sched_getaffinity(0)):
        p.error("Worker count exceeds this job's CPU affinity; request more CPUs")
    log(f"START workers={args.workers}, affinity_cpus={len(os.sched_getaffinity(0))}, DESI={base.DESISPEC_VERSION}")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    manifest, a = prepare_selection(args)
    if args.prepare_only:
        return
    signature = dict(manifest=str(manifest), manifest_size=manifest.stat().st_size,
                     manifest_mtime_ns=manifest.stat().st_mtime_ns,
                     rows=len(a), workers=args.workers, desi=base.DESISPEC_VERSION,
                     filter=base.FILTER_DESCRIPTION, chunk_rows=CHUNK, output=str(args.output),
                     script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.signature = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
    mode = f"smoke_{args.smoke_groups}" if args.smoke_groups else "full"
    directory = args.work_dir / mode
    directory.mkdir(exist_ok=True)
    config_path = directory / "build_config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != signature:
        raise RuntimeError("Resume configuration mismatch; do not mix builds")
    atomic_json(config_path, signature)
    blocks = (len(a) + CHUNK - 1) // CHUNK
    n_shards = min(args.workers, blocks)
    boundaries = [min((i * blocks // n_shards) * CHUNK, len(a)) for i in range(n_shards)] + [len(a)]
    tasks = [(str(manifest), boundaries[i], boundaries[i + 1], str(directory / f"shard_{i:03d}.h5"),
              str(args.zall), str(args.coadd_root)) for i in range(n_shards)]
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=CTX) as pool:
        for result in pool.map(build_shard, tasks):
            log(f"shard complete: {result}")
    output = merge_shards(args, tasks, len(a), len(group_slices(a)))
    log("Starting parallel full chunk and metadata verification")
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=CTX) as pool:
        checks = list(pool.map(validate_shard, [(t, str(output)) for t in tasks]))
    assert sum(c["rows"] for c in checks) == len(a)
    if len(np.unique(a["TARGETID"])) != len(a):
        raise RuntimeError("Duplicate output TARGETID")
    with h5py.File(output, "r+") as h:
        assert all(h[k].shape[0] == len(a) for k in base.ZALL_FIELDS + tuple(base.SPECTRUM_FIELDS))
        h.attrs["VALIDATION_COMPLETE"] = True
        h.attrs["BUILD_COMPLETE"] = True
        h.flush()
    report = dict(output=str(args.output), rows=len(a), unique_targetids=len(a), workers=args.workers,
                  filter=base.FILTER_DESCRIPTION, desispec_version=base.DESISPEC_VERSION,
                  validation="all metadata and all compressed spectral chunks exact; source endpoint recomputation",
                  shards=checks)
    atomic_json(directory / "validation.json", report)
    if output != args.output:
        if args.output.exists():
            raise FileExistsError(args.output)
        output.replace(args.output)
    log(f"COMPLETE {args.output}: {len(a):,} unique spectra; validation passed")


if __name__ == "__main__":
    main()
