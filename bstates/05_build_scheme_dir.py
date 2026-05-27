#!/usr/bin/env python3
"""
05_build_scheme_dir.py
============================================================================
Assemble a WESTPA basis-state directory for a chosen seeding scheme.

Given a list of ORIGINAL struct indices (the struct_X numbering used in
03_bstates_final/ and in 04_analyze) and a user-defined output directory name,
this tool:

  1. copies each selected 03_bstates_final/StructureFiles/struct_X directory
     into  <OUTDIR>/StructureFiles/ ,
  2. RENAMES the copies sequentially from struct_0 (first selected -> struct_0,
     second -> struct_1, ...),  preserving the order the indices are given,
  3. writes <OUTDIR>/bstates.txt  re-indexed 0..N-1 with a uniform weight 1/N
     and the canonical  StructureFiles/struct_i  path (the same layout WESTPA
     reads in 03_bstates_final/bstates.txt),
  4. writes <OUTDIR>/struct_mapping.log recording  new_index <- orig struct_X
     <- bstate_id  so every renamed seed is traceable to its source.

Layout produced (canonical, matches 03_bstates_final/):

    <OUTDIR>/
        bstates.txt
        struct_mapping.log
        StructureFiles/
            struct_0/   (copy of orig struct_<first>)
            struct_1/   (copy of orig struct_<second>)
            ...

USAGE
    # explicit indices, in the order you want them renamed
    python3 05_build_scheme_dir.py --name Scheme_E \
        --structs 4 11 35 50 64 52 60 65 42 53 37 22 31

    # or pull the indices straight from a recommendations CSV (one scheme)
    python3 05_build_scheme_dir.py --name Scheme_E \
        --from-csv 04_analyze/results/recommendations_schemeE.csv

Options:
    --name STR          output directory name (created under repo root unless
                        an absolute/relative path is given)
    --structs N [N...]  original struct indices, in desired order
    --from-csv FILE     read struct_index column from a recommendations CSV
                        (use --scheme to filter if the CSV has many schemes)
    --scheme STR        when using --from-csv, keep only rows whose 'scheme'
                        column == STR  (default: all rows)
    --source DIR        source bstates dir (default: 03_bstates_final)
    --force             overwrite OUTDIR if it already exists
============================================================================
"""
import argparse
import csv
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load_mapping(source):
    """orig struct_index -> bstate_id, from <source>/struct_mapping.log."""
    path = os.path.join(source, "struct_mapping.log")
    mapping = {}
    if not os.path.isfile(path):
        return mapping
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                sidx = int(parts[0].replace("struct_", ""))
            except (ValueError, IndexError):
                continue
            mapping[sidx] = parts[1] if len(parts) > 1 else "?"
    return mapping


def structs_from_csv(path, scheme):
    """Read struct indices (preserving file order) from a recommendations CSV."""
    out = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if "struct_index" not in row:
                sys.exit(f"ERROR: {path} has no 'struct_index' column.")
            if scheme and row.get("scheme", "") != scheme:
                continue
            out.append(int(row["struct_index"]))
    return out


def main():
    ap = argparse.ArgumentParser(description="Build a WESTPA bstate directory for a seeding scheme.")
    ap.add_argument("--name", required=True, help="output directory name (or path)")
    ap.add_argument("--structs", type=int, nargs="+", help="original struct indices, in desired order")
    ap.add_argument("--from-csv", help="recommendations CSV to read struct_index from")
    ap.add_argument("--scheme", default=None, help="filter CSV rows to this scheme value")
    ap.add_argument("--source", default=os.path.join(HERE, "03_bstates_final"),
                    help="source bstates dir (default: 03_bstates_final)")
    ap.add_argument("--force", action="store_true", help="overwrite OUTDIR if it exists")
    args = ap.parse_args()

    # ---- resolve the ordered list of source struct indices ----------------
    if args.structs and args.from_csv:
        sys.exit("ERROR: give either --structs or --from-csv, not both.")
    if args.structs:
        structs = list(args.structs)
    elif args.from_csv:
        structs = structs_from_csv(args.from_csv, args.scheme)
    else:
        sys.exit("ERROR: provide --structs N [N...] or --from-csv FILE.")
    if not structs:
        sys.exit("ERROR: no struct indices selected.")

    # de-duplicate while preserving order (a struct seeded twice is a no-op)
    seen, ordered = set(), []
    for s in structs:
        if s in seen:
            print(f"  WARNING: struct_{s} listed more than once; keeping first occurrence.")
            continue
        seen.add(s)
        ordered.append(s)
    structs = ordered

    source = os.path.abspath(args.source)
    src_sf = os.path.join(source, "StructureFiles")
    if not os.path.isdir(src_sf):
        sys.exit(f"ERROR: source StructureFiles dir not found: {src_sf}")

    # ---- pre-flight: every requested struct must exist --------------------
    missing = [s for s in structs if not os.path.isdir(os.path.join(src_sf, f"struct_{s}"))]
    if missing:
        sys.exit(f"ERROR: these struct dirs are missing in {src_sf}: "
                 + ", ".join(f"struct_{s}" for s in missing))

    bstate_of = load_mapping(source)

    # ---- prepare output dir -----------------------------------------------
    outdir = args.name if os.path.isabs(args.name) else os.path.join(HERE, args.name)
    if os.path.exists(outdir):
        if not args.force:
            sys.exit(f"ERROR: {outdir} already exists. Re-run with --force to overwrite.")
        print(f"  --force: removing existing {outdir}")
        shutil.rmtree(outdir)
    out_sf = os.path.join(outdir, "StructureFiles")
    os.makedirs(out_sf)

    # ---- copy + rename sequentially from 0 --------------------------------
    n = len(structs)
    weight = 1.0 / n
    rows = []   # (new_index, orig_index, bstate_id)
    print(f"Building '{outdir}'  ({n} seeds, uniform weight {weight:.12e})")
    for new_i, orig in enumerate(structs):
        src = os.path.join(src_sf, f"struct_{orig}")
        dst = os.path.join(out_sf, f"struct_{new_i}")
        shutil.copytree(src, dst)
        bid = bstate_of.get(orig, "?")
        rows.append((new_i, orig, bid))
        print(f"  struct_{new_i:<3d} <- struct_{orig:<3d}  ({bid})")

    # ---- bstates.txt (canonical StructureFiles/struct_i layout) -----------
    bpath = os.path.join(outdir, "bstates.txt")
    with open(bpath, "w") as fh:
        for new_i, _, _ in rows:
            fh.write(f"{new_i}    {weight:.17e}    StructureFiles/struct_{new_i}\n")

    # ---- struct_mapping.log -----------------------------------------------
    mpath = os.path.join(outdir, "struct_mapping.log")
    with open(mpath, "w") as fh:
        fh.write("# ===========================================================================\n")
        fh.write(f"# Scheme bstate directory: {os.path.basename(outdir.rstrip('/'))}\n")
        fh.write(f"# Source                 : {source}\n")
        fh.write(f"# Seeds                  : {n}  (uniform weight {weight:.12e})\n")
        fh.write("# Each seed was copied from the source StructureFiles/ and renamed\n")
        fh.write("# sequentially from struct_0. Columns map the new index back to the\n")
        fh.write("# original struct_X numbering (03_bstates_final / 04_analyze) and bstate_id.\n")
        fh.write("#\n")
        fh.write("# new_index   new_dir        orig_struct      bstate_id\n")
        fh.write("# ===========================================================================\n")
        for new_i, orig, bid in rows:
            fh.write(f"{new_i:<11d} struct_{new_i:<8d} struct_{orig:<10d} {bid}\n")

    print(f"\nWrote:\n  {bpath}\n  {mpath}\n  {out_sf}/  ({n} struct_* dirs)")


if __name__ == "__main__":
    main()
