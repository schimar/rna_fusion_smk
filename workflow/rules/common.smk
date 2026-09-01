import glob
from pathlib import Path
import os
import shutil
import re
import pandas as pd
from snakemake.utils import validate

# -----------------------------------------------------------------------------
# Configuration variables
# -----------------------------------------------------------------------------
runid = config['runid']
bcldir = config['bcldir']
analysis_path = config['analysis_path']
rgid = bcldir.split('_')[1]  # config['rgID']
bed_file = config['bed']
bed_file_stem = bed_file.split('.')[0]

# Reference genome (configurable via config.yaml 'ref' = file stem):
#   resources/{ref}.fasta | .gtf | .txt (chrom sizes) | .dict | .fasta.fai
#   star index dir: resources/star_{ref}/
ref = config['ref']
genome_fa = f"resources/{ref}.fasta"
genome_gtf = f"resources/{ref}.gtf"
genome_txt = f"resources/{ref}.txt"
genome_dict = f"resources/{ref}.dict"
genome_fai = f"resources/{ref}.fasta.fai"
star_idx = f"resources/star_{ref}"
#r1_read_structure = config["r1_read_structure"]
#r2_read_structure = config["r2_read_structure"]
samples_tsv = '/'.join([runid, 'samples.tsv'])
sample_sheet = bcldir + '/SampleSheet_rna.csv'

# Create necessary directories
pr = Path('/'.join([runid, 'results/']))
pr.mkdir(parents=True, exist_ok=True)
pl = Path('/'.join([runid, 'logs/']))
pl.mkdir(parents=True, exist_ok=True)


def parse_sample_sheet(sheet_path):
    """
    Parse a SampleSheet CSV (v1 or v2) and return the ordered list of
    Sample_ID values.

    v1: header row has 'Sample_ID' at column 0.
    v2 (e.g. NovaSeqX+): 'Sample_ID' may be in any column; the header row is
    located by scanning for the row containing 'Sample_ID'.
    """
    sample_ids = []
    if not os.path.isfile(sheet_path):
        print(f"Warning: SampleSheet not found: {sheet_path}")
        return sample_ids
    col = None
    with open(sheet_path) as fh:
        for line in fh:
            line = line.strip('\r\n')
            if not line:
                continue
            cols = [c.strip() for c in line.split(',')]
            if col is None:
                # header row: the row that contains the 'Sample_ID' column
                if 'Sample_ID' in cols:
                    col = cols.index('Sample_ID')
                continue
            # data rows
            if col < len(cols) and cols[col]:
                sample_ids.append(cols[col])
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
Path('/'.join([runid])).mkdir(parents=True, exist_ok=True)
samples_df.to_csv(samples_tsv, sep='\t', index=False)
print(f"WTS samples: {wts_samples}")
print(f"sample table written to {samples_tsv}")

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
