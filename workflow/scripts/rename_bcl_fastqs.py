#!/usr/bin/env python3
"""Rename bcl-convert fastq outputs.

bcl-convert emits files named  <sample>_S<index>_R<read>_001.fastq.gz
(trailing _S<index>_R<read>_001.fastq.gz is appended by bcl-convert).

This script moves ALL demuxed fastqs (including non-WTS samples) out of the
fresh tmp dir into the final fastq dir, renaming them to
  <sample>.R<read>_001.fastq.gz
so that a later `rm -rf tmp` does not destroy demultiplexed data that other
(parts of future) workflows may still need.

Usage:
    python3 rename_bcl_fastqs.py <bcl_tmp_dir> <outdir>
"""

import glob
import os
import re
import sys

PATTERN = re.compile(r"^(.*)_S[0-9]+_R([12])_001\.fastq\.gz$")


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: rename_bcl_fastqs.py <tmp_dir> <outdir>", file=sys.stderr)
        return 2
    bcl_tmp, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    fastqs = sorted(glob.glob(os.path.join(bcl_tmp, "*.fastq.gz")))
    moves = []
    errors = []
    for fq in fastqs:
        base = os.path.basename(fq)
        m = PATTERN.match(base)
        if not m:
            errors.append(f"unexpected fastq name {base!r}")
            continue
        dest = os.path.join(outdir, f"{m.group(1)}.R{m.group(2)}_001.fastq.gz")
        if os.path.exists(dest):
            errors.append(f"destination already exists: {dest}")
            continue
        moves.append((fq, dest))

    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        print("no fastq files were moved; preserving bcl-convert output", file=sys.stderr)
        return 1

    for fq, dest in moves:
        os.replace(fq, dest)
    print(f"moved {len(moves)} fastq files into {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())