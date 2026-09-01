#!/usr/bin/env python3
"""GGTF -> flat BED6 of rRNA + Mt_rRNA loci (SPEC §12.3)."""
import argparse, sys, pandas as pd

RRNA_BIOTYPES = {"rrna", "mt_rrna"}

def parse_gtf_attr(s):
    out = {}
    for part in s.rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        k, _, v = part.partition(" ")
        out.setdefault(k.strip(), []).append(v.strip().strip('"'))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--gtf", required=True)
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    rows = []
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            a = parse_gtf_attr(c[8])
            biotypes = {b.lower() for b in a.get("gene_biotype", [])}
            if not (biotypes & RRNA_BIOTYPES):
                continue
            name = (a.get("gene_name") or ["rRNA"])
            name = name[0] if isinstance(name, list) else name
            rows.append((c[0], max(0, int(c[3]) - 1), int(c[4]), name, c[6] if len(c) > 6 else "."))

    if not rows:
        print("WARNING: no rRNA features found in " + args.gtf +
              " (gene_biotype in {rRNA, Mt_rRNA}); BED is empty — check GTF "
              "completeness (rDNA repeat units often absent/collapsed).", file=sys.stderr)

    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    merged = []
    for r in rows:
        if merged and merged[-1][0] == r[0] and r[1] <= merged[-1][2]:
            merged[-1] = (merged[-1][0], merged[-1][1], max(merged[-1][2], r[2]), merged[-1][3], merged[-1][4])
        else:
            merged.append(r)

    with open(args.out, "w")as out:
        for chrom, s, e, name, strand in merged:
            out.write(f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}\n")
    print(f"Wrote {len(merged)} rRNA interval(s) to {args.out}")

if __name__ == "__main__":
    main()
