#!/usr/bin/env python3
"""Create a MultiQC custom-content table from per-sample rRNA QC results."""
import csv
import json


data = {}
for path in snakemake.input:
    with open(path, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"), None)
    if row is None:
        continue
    data[row["sample"]] = {
        "rRNA fraction (%)": float(row["rrna_fraction_percent"]),
        "rRNA reads": int(row["rrna_reads"]),
        "Primary mapped reads": int(row["mapped_primary_reads"]),
    }

content = {
    "id": "rrna_contamination",
    "section_name": "rRNA contamination",
    "description": "Primary mapped alignments overlapping rRNA loci.",
    "plot_type": "table",
    "pconfig": {
        "id": "rrna_contamination_table",
        "headers": {
            "rRNA fraction (%)": {"format": "{:.2f}"},
            "rRNA reads": {"format": "{:,d}"},
            "Primary mapped reads": {"format": "{:,d}"},
        },
    },
    "data": data,
}

with open(snakemake.output[0], "w") as handle:
    json.dump(content, handle, indent=2)
    handle.write("\n")