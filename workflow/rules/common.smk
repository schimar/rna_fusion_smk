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
no_lane_splitting = str(not config.get('lane_splitting', False)).lower()


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
        return []
    blocks = []
    current_block = None
    with open(sheet_path, newline="") as fh:
        for row in csv.reader(fh):
            cols = [cell.strip() for cell in row]
            if not any(cols):
                continue
            if 'Sample_ID' in cols:
                current_block = {
                    "column": cols.index('Sample_ID'),
                    "sample_ids": [],
                }
                blocks.append(current_block)
                continue

            if current_block is None:
                continue
            column = current_block["column"]
            if column < len(cols) and cols[column]:
                current_block["sample_ids"].append(cols[column])

    if not blocks:
        print(f"Warning: Sample_ID header not found in {sheet_path}")
        return []

    selected_block = max(
        blocks,
        key=lambda block: sum('WTS' in sample_id for sample_id in block["sample_ids"]),
    )
    sample_ids = selected_block["sample_ids"]
    print(f"Parsed {len(sample_ids)} sample(s) from {sheet_path}")
    return sample_ids


# -----------------------------------------------------------------------------
# Sample definition & selection
# The SampleSheet is the source of truth (parsed at parse time). A sample is a
# WTS sample iff its full Sample_ID contains the substring 'WTS'. No other part
# of the ID is interpreted (kept whole, including any _S<n> suffix).
# -----------------------------------------------------------------------------
all_sample_ids = parse_sample_sheet(sample_sheet)
wts_samples = [sid for sid in all_sample_ids if 'WTS' in sid]

if not wts_samples:
    print("WARNING: no WTS samples found in SampleSheet " + sample_sheet)

idkeys = list(wts_samples)

# Write the derived sample table (all samples + WTS flag) as a transparency
# artifact. This is derived, NOT the source of truth.
samples_df = pd.DataFrame({'sample_id': all_sample_ids})
samples_df['is_wts'] = samples_df['sample_id'].apply(lambda s: 'WTS' in s)
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

# -----------------------------------------------------------------------------
# Workflow hooks
# -----------------------------------------------------------------------------
onstart:
    shell("gitrev=$(git rev-parse HEAD) && echo \"--------------------------------------------- \n Running rna_fusion workflow, with \n git branch $gitrev \n log file in \n {log} \n & run data in \n {runid} \n  ---------------------------------------------\" 2>&1 | tee {log}")
onsuccess:
    shell("cp -v {log} {runid}/logs/")
onsuccess:
    shutil.rmtree(".snakemake")
onerror:
    shell("cp -v {log} {runid}/logs/")
