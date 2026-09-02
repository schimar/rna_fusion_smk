#!/usr/bin/env python3
"""Rewrite GTF sequence names to the matching reference FASTA contig IDs."""
import argparse
import sys


def fasta_contigs(path):
    contigs = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                contigs.add(line[1:].split(None, 1)[0])
    return contigs


def matching_contig(contig, reference_contigs):
    for candidate in (contig, f"chr{contig}", "chrM" if contig == "MT" else ""):
        if candidate in reference_contigs:
            return candidate
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-o", "--out", required=True)
    args = parser.parse_args()

    reference_contigs = fasta_contigs(args.fasta)
    written = 0
    skipped = 0
    with open(args.gtf) as source, open(args.out, "w") as destination:
        for line in source:
            if line.startswith("#"):
                destination.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                destination.write(line)
                continue
            contig = matching_contig(fields[0], reference_contigs)
            if contig is None:
                skipped += 1
                continue
            fields[0] = contig
            destination.write("\t".join(fields) + "\n")
            written += 1

    if not written:
        raise RuntimeError("No GTF records matched reference FASTA contigs")
    print(f"Wrote {written} records; skipped {skipped} unmatched records", file=sys.stderr)


if __name__ == "__main__":
    main()