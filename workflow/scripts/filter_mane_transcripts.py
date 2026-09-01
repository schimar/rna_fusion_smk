#!/usr/bin/env python3
"""Write genome-wide MANE_Select and MANE_Plus_Clinical transcript records."""
import argparse
import sys


MANE_TAGS = {"MANE_Select", "MANE_Plus_Clinical"}


def parse_attrs(attributes):
    parsed = {}
    for part in attributes.rstrip(";").split(";"):
        key, _, value = part.strip().partition(" ")
        if key:
            parsed.setdefault(key, []).append(value.strip().strip('"'))
    return parsed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True)
    parser.add_argument("-o", "--out", required=True)
    args = parser.parse_args()

    transcript_ids = set()
    with open(args.gtf) as source:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attributes = parse_attrs(fields[8])
            tags = {
                tag
                for value in attributes.get("tag", [])
                for tag in value.split(",")
            }
            if tags & MANE_TAGS:
                transcript_ids.update(attributes.get("transcript_id", []))

    if not transcript_ids:
        print(
            "WARNING: no MANE_Select or MANE_Plus_Clinical transcripts found in "
            + args.gtf,
            file=sys.stderr,
        )

    with open(args.gtf) as source, open(args.out, "w") as destination:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attributes = parse_attrs(fields[8])
            if set(attributes.get("transcript_id", [])) & transcript_ids:
                destination.write(line)

    print(f"Selected {len(transcript_ids)} MANE transcript(s)")


if __name__ == "__main__":
    main()