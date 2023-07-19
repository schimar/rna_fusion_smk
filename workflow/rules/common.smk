import glob
from pathlib import Path

import pandas as pd
#from snakemake.remote import FTP
from snakemake.utils import validate

#ftp = FTP.RemoteProvider()
fqDir = ""   #"/mnt/illumina/230209_A01272_0035_BHTVFGDRX2"
#fqDir = "/mnt/illumina/development/bench_data/rna/fq"
runid = config['runID']
rgid =  config['rgID']

validate(config, schema="../../config/schemas/config.schema.yaml")


units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../../config/schemas/units.schema.yaml")

samples_dict = dict(zip(units['sample_name'], zip(units['fq'])))  #, units['fq2'])))



pth = Path(units['fq'][0])

fqDir = str(pth.parent)

idkeys = list(samples_dict.keys())

# NOTE: add a quick info about regex below!!!
common_constraint = "[0-9A-Za-z\-\_]+[^\/L][^\/umi][^\/chr7][^\/SE]"

# read structure of our UMIs 
r1_read_structure = "5M2S+T"  #"8M+T"
r2_read_structure = "5M2S+T"  #"8M+T"




