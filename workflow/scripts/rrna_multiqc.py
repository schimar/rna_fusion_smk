#!/usr/bin/env python3
"""Create a MultiQC custom-content General Stats block from per-sample rRNA QC results."""
import csv
import json


data = {}
for path in snakemake.input:
    with open(path, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"), None)
    if row is None:
        continue
    data[row["sample"]] = {
        "rrna_fraction_pct": float(row["rrna_fraction_percent"]),
        "rrna_reads": int(row["rrna_reads"]),
        "mapped_primary_reads": int(row["mapped_primary_reads"]),
    }

content = {
    "id": "rrna_contamination",
    "section_name": "rRNA contamination",
    "description": "Primary mapped alignments overlapping rRNA loci.",
    "plot_type": "generalstats",
    "headers": {
        "rrna_fraction_pct": {"title": "rRNA fraction (%)", "format": "{:.2f}"},
        "rrna_reads": {"title": "rRNA reads", "format": "{:,d}", "hidden": True},
        "mapped_primary_reads": {
            "title": "Primary mapped reads",
            "format": "{:,d}",
            "hidden": True,
        },
    },
    "data": data,
}

with open(snakemake.output[0], "w") as handle:
    json.dump(content, handle, indent=2)
    handle.write("\n")