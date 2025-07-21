import glob
from pathlib import Path
import os
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
r1_read_structure = config["r1_read_structure"]
r2_read_structure = config["r2_read_structure"]
units_wd = '/'.join([runid, config["units"]])

# Create necessary directories
pr = Path('/'.join([runid, 'results/']))
pr.mkdir(parents=True, exist_ok=True)
pl = Path('/'.join([runid, 'logs/']))
pl.mkdir(parents=True, exist_ok=True)

def detect_lanes_and_create_units(bcl2fq_path, output_path, file_pattern_suffix="R1_001.fastq.gz"):
    """
    Detect available lanes in the dataset and create a units.tsv file
    
    Parameters:
    -----------
    bcl2fq_path : str
        Path to the bcl2fq output directory
    output_path : str
        Path where units.tsv will be saved
    file_pattern_suffix : str
        Suffix pattern to use for finding R1 files
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing sample information
    """
    lsda = os.listdir(bcl2fq_path)
    
    # Find all R1 files (excluding Undetermined)
    r1_files = [f for f in lsda if f.endswith(file_pattern_suffix) and 'Undetermined' not in f]
    
    if not r1_files:
        print(f"Warning: No R1 files found with pattern '*{file_pattern_suffix}' in {bcl2fq_path}")
        return pd.DataFrame(columns=['sample_name', 'unit_name', 'fq'])
    
    # Check which lanes are present in the dataset by examining all files
    available_lanes = set()
    for file in r1_files:
        # Look for lane patterns L001, L002, L003, L004
        lane_matches = re.findall('L00[1-4]', file)
        if lane_matches:
            available_lanes.update(lane_matches)
    
    # If no lanes detected or ambiguous, assume L001 only
    if not available_lanes:
        print("No lane information found in filenames, defaulting to L001")
        lanes = ['L001' for _ in r1_files]
    else:
        # Extract lanes for each file
        lanes = []
        for file in r1_files:
            lane_match = re.findall('L00[1-4]', file)
            lanes.append(lane_match[0] if lane_match else 'L001')
    
    # Extract sample names from filenames
    sample_names = ['_'.join(x.split('_')[0:2]) for x in r1_files]
    
    # Create full file paths
    fq_paths = ['/'.join([runid, 'results/bcl2fq', fqgz]) for fqgz in r1_files]
    
    # Create DataFrame
    df = pd.DataFrame({'sample_name': sample_names, 'unit_name': lanes, 'fq': fq_paths})
    df.sort_values(['sample_name', 'unit_name'], ascending=[True, True], inplace=True)
    
    # Print information about detected lanes
    detected_lanes = sorted(available_lanes)
    print(f"Detected lanes: {', '.join(detected_lanes)}")
    
    # Save to file
    df.to_csv(output_path, sep='\t', index=False)
    
    return df

# -----------------------------------------------------------------------------
# Load or create units file
# -----------------------------------------------------------------------------
units_file_exists = Path(units_wd).is_file()
units_file_has_data = False

if units_file_exists:
    # Check if the file has data beyond headers
    try:
        units_df = pd.read_csv(units_wd, sep="\t")
        units_file_has_data = len(units_df) > 0
    except:
        units_file_has_data = False

if units_file_exists and units_file_has_data:
    # File exists and has data
    units = (pd.read_csv(units_wd, sep="\t", dtype={"sample_name": str, "unit_name": str})
        .set_index(["sample_name", "unit_name"], drop=False)
        .sort_index()
    )
    validate(units, schema="../../config/schemas/units.schema.yaml")
else:
    # File doesn't exist, is empty, or has only headers
    # Create bcl2fq results directory if needed
    pr = Path('/'.join([runid, 'results/bcl2fq/']))
    pr.mkdir(parents=True, exist_ok=True)
    
    # Determine the appropriate file pattern suffix based on what's in the directory
    bcl2fq_dir = '/'.join([runid, "results/bcl2fq/"])
    
    # First try with standard Illumina naming
    file_pattern = "R1_001.fastq.gz"
    
    # If no files match the standard pattern, try the FDA challenge pattern
    if not any(f.endswith(file_pattern) for f in os.listdir(bcl2fq_dir) if 'Undetermined' not in f):
        file_pattern = "R1.fq.gz"
    
    # Detect lanes and create units file
    units = detect_lanes_and_create_units(
        bcl2fq_path=bcl2fq_dir,
        output_path='/'.join([runid, 'units.tsv']),
        file_pattern_suffix=file_pattern
    )
    
    validate(units, schema="../../config/schemas/units.schema.yaml")

# Create samples dictionary and ID keys
samples_dict = dict(zip(units['sample_name'], zip(units['fq'])))
idkeys = list(samples_dict.keys())

# -----------------------------------------------------------------------------
# Wildcard constraints
# -----------------------------------------------------------------------------
common_constraint = r"(?!.*dedupe)[^/]+"

# -----------------------------------------------------------------------------
# Workflow hooks
# -----------------------------------------------------------------------------
onstart:
    shell("gitrev=$(git rev-parse HEAD) && echo \"--------------------------------------------- \n Running rna_fusion workflow, with \n git branch $gitrev \n log file in \n {log} \n & run data in \n {runid} \n  ---------------------------------------------\" 2>&1 | tee {log}")
onsuccess:
    shell("cp -v {log} {runid}/logs/")
onerror:
    shell("cp -v {log} {runid}/logs/")
