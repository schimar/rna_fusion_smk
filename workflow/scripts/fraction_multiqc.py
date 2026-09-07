#!/usr/bin/env python3
"""Create a MultiQC custom-content General Stats block from per-sample fraction QC results.

Parameterized via snakemake.params:
    metric: short id used for the JSON "id" and the data key (e.g. "rrna", "globin", "mane")
    column: fraction column name in the per-sample TSV (e.g. "rrna_fraction_percent")
    title:  display title for the generalstats column
    description: description for the generalstats column
"""
import csv
import json


metric = snakemake.params["metric"]
column = snakemake.params["column"]
title = snakemake.params["title"]
description = snakemake.params["description"]

data = {}
for path in snakemake.input:
    with open(path, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"), None)
    if row is None:
        continue
    data[row["sample"]] = {metric: float(row[column])}

content = {
    "id": f"{metric}_fraction",
    "section_name": title,
    "description": description,
    "plot_type": "generalstats",
    "headers": {
        metric: {
            "title": title,
            "description": description,
            "format": "{:.2f}",
            "suffix": "%",
        },
    },
    "data": data,
}

with open(snakemake.output[0], "w") as handle:
    json.dump(content, handle, indent=2)
    handle.write("\n")