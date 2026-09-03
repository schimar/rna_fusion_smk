#!/usr/bin/env python3
"""Rename bcl-convert fastq outputs.

bcl-convert emits files named either:
    <sample>_S<index>_R<read>_001.fastq.gz
    <sample>_S<index>_L<lane>_R<read>_001.fastq.gz

This script moves ALL demuxed fastqs (including non-WTS samples) out of the
fresh tmp dir into the final fastq dir, renaming them to either
    <sample>.R<read>_001.fastq.gz
    <sample>_L<lane>.R<read>_001.fastq.gz
so that a later `rm -rf tmp` does not destroy demultiplexed data that other
(parts of future) workflows may still need.

Per-file semantics (not all-or-nothing):
  - A file that matches PATTERN is moved, overwriting any stale destination
    from a previous run of this rule. A rerun with an updated sample sheet
    must not be blocked by leftover output from an earlier run.
  - A file that does NOT match PATTERN is left in place (not moved) and
    reported as an error. Any such mismatch fails the whole script (non-zero
    exit) so a naming assumption violation is never silently swallowed -- but
    it no longer prevents every OTHER, correctly-named file from being moved
    first. Fail loud, don't fail wide.

Usage:
    python3 rename_bcl_fastqs.py <bcl_tmp_dir> <outdir>
"""
import glob
import os
import re
import sys

PATTERN = re.compile(r"^(.*)_S[0-9]+(_L[0-9]+)?_R([12])_001\.fastq\.gz$")


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: rename_bcl_fastqs.py <tmp_dir> <outdir>", file=sys.stderr)
        return 2

    bcl_tmp, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    fastqs = sorted(glob.glob(os.path.join(bcl_tmp, "*.fastq.gz")))
    moved = []
    errors = []

    for fq in fastqs:
        base = os.path.basename(fq)
        m = PATTERN.match(base)
        if not m:
            errors.append(f"unexpected fastq name {base!r} (left in place, not moved)")
            continue

        lane = m.group(2) or ""
        dest = os.path.join(outdir, f"{m.group(1)}{lane}.R{m.group(3)}_001.fastq.gz")

        # Overwrite any stale destination from a previous run. A pre-existing
        # file here is expected on rerun (e.g. sample sheet gained new rows)
        # and must not block moving files that are otherwise fine.
        try:
            os.replace(fq, dest)
            moved.append((fq, dest))
        except OSError as exc:
            errors.append(f"failed to move {base!r} -> {dest!r}: {exc}")

    print(f"moved {len(moved)} fastq files into {outdir}")

    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        print(
            f"{len(errors)} file(s) could not be moved -- see errors above",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

