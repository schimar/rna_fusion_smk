import glob
import csv
from pathlib import Path
import os
import shutil
import re
import pandas as pd
from snakemake.utils import validate

# -----------------------------------------------------------------------------
# Configuration variables
# -----------------------------------------------------------------------------
bcldir = config['bcldir']
runid = config.get('runid') or Path(bcldir).name
analysis_path = config['analysis_path']
final_dest = config.get('final_dest') or (
    "/home/schilling_m1/smb/Analyses/00_Tests/vc_rna/" + Path(bcldir).name
)
rgid = config.get('rgid') or Path(bcldir).name
bed_file = config['bed']
bed_file_stem = bed_file.split('.')[0]

# Reference genome (configurable via config.yaml 'ref' = file stem):
#   resources/{ref}.fasta | .gtf | .txt (chrom sizes) | .dict | .fasta.fai
#   star index dir: resources/star_{ref}/
ref = config['ref']
genome_fa = f"resources/{ref}.fasta"
genome_gtf_source = f"resources/{ref}.gtf"
genome_gtf = f"resources/{ref}.normalized.gtf"
genome_txt = f"resources/{ref}.txt"
genome_dict = f"resources/{ref}.dict"
genome_fai = f"resources/{ref}.fasta.fai"
star_idx = f"resources/star_{ref}"
#r1_read_structure = config["r1_read_structure"]
#r2_read_structure = config["r2_read_structure"]
samples_tsv = '/'.join([runid, 'samples.tsv'])
sample_sheet_config = config.get('sample_sheet') or 'SampleSheet.csv'
sample_sheet = (
    sample_sheet_config
    if os.path.isabs(sample_sheet_config)
    else os.path.join(bcldir, sample_sheet_config)
)
lane_splitting = config.get('lane_splitting', False)
no_lane_splitting = str(not lane_splitting).lower()


def parse_sample_sheet(sheet_path):
    """
    Parse a SampleSheet CSV (v1 or v2) and return the ordered list of
    Sample_ID values.

    v1: header row has 'Sample_ID' at column 0.
    v2 (e.g. NovaSeqX+): 'Sample_ID' may be in any column; the header row is
    located by scanning for the row containing 'Sample_ID'.
    """
    if not os.path.isfile(sheet_path):
        print(f"Warning: SampleSheet not found: {sheet_path}")
        return [], []
    blocks = []
    current_block = None
    with open(sheet_path, newline="") as fh:
        for row in csv.reader(fh):
            cols = [cell.strip() for cell in row]
            if not any(cols):
                continue
            if len(cols) == 1 and cols[0].startswith("[") and cols[0].endswith("]"):
                current_block = None
                continue
            if 'Sample_ID' in cols:
                current_block = {
                    "header": cols,
                    "rows": [],
                }
                blocks.append(current_block)
                continue

            if current_block is None:
                continue
            if len(cols) >= len(current_block["header"]):
                current_block["rows"].append(cols)

    if not blocks:
        print(f"Warning: Sample_ID header not found in {sheet_path}")
        return [], []

    selected_block = max(
        blocks,
        key=lambda block: sum(
            'WTS' in row[block["header"].index('Sample_ID')]
            for row in block["rows"]
        ),
    )
    header = selected_block["header"]
    sample_column = header.index('Sample_ID')
    samples = []
    for row in selected_block["rows"]:
        sample_id = row[sample_column]
        if sample_id:
            samples.append((sample_id, row))
    print(f"Parsed {len(samples)} sample(s) from {sheet_path}")
    return samples, header


# -----------------------------------------------------------------------------
# Sample definition & selection
# The SampleSheet is the source of truth (parsed at parse time). A sample is a
# WTS sample iff its full Sample_ID contains the substring 'WTS'. No other part
# of the ID is interpreted (kept whole, including any _S<n> suffix).
# -----------------------------------------------------------------------------
sample_rows, sample_header = parse_sample_sheet(sample_sheet)
all_sample_ids = [sample_id for sample_id, _ in sample_rows]

if lane_splitting:
    if 'Lane' not in sample_header:
        raise ValueError("lane_splitting=true requires a Lane column in the SampleSheet")
    lane_column = sample_header.index('Lane')
    sample_units = []
    for sample_id, row in sample_rows:
        lane_values = row[lane_column].split('+') if row[lane_column] else []
        if not lane_values:
            raise ValueError(f"Missing Lane value for Sample_ID {sample_id!r}")
        for lane_value in lane_values:
            if not lane_value.isdigit():
                raise ValueError(f"Invalid Lane value {lane_value!r} for Sample_ID {sample_id!r}")
            sample_units.append((sample_id, f"L{int(lane_value):03d}"))
else:
    sample_units = [(sample_id, "") for sample_id, _ in sample_rows]

wts_samples = list(dict.fromkeys(
    f"{sample_id}_{lane}" if lane else sample_id
    for sample_id, lane in sample_units
    if 'WTS' in sample_id
))

if not wts_samples:
    print("WARNING: no WTS samples found in SampleSheet " + sample_sheet)

idkeys = list(wts_samples)

# Write the derived sample table (all samples + WTS flag) as a transparency
# artifact. This is derived, NOT the source of truth.
samples_df = pd.DataFrame(sample_units, columns=['sample_id', 'lane'])
samples_df['analysis_unit'] = samples_df.apply(
    lambda row: f"{row.sample_id}_{row.lane}" if row.lane else row.sample_id,
    axis=1,
)
samples_df['is_wts'] = samples_df['sample_id'].apply(lambda sample_id: 'WTS' in sample_id)
print(f"WTS samples: {wts_samples}")


rule write_samples_table:
    output:
        samples_tsv,
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        samples_df.to_csv(output[0], sep='\t', index=False)

# -----------------------------------------------------------------------------
# Wildcard constraints
# -----------------------------------------------------------------------------
common_constraint = r"(?!.*dedupe)[^/]+"    #(?!.*cons)
