import glob
from pathlib import Path
import os
import re
import pandas as pd
#from snakemake.remote import FTP
from snakemake.utils import validate

#ftp = FTP.RemoteProvider()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#fqDir = ""   #"/mnt/illumina/230209_A01272_0035_BHTVFGDRX2"
runid = config['runid']
bcldir = config['bcldir']
rgid =  config['rgID']
# read structure of our UMIs 
r1_read_structure = config["r1_read_structure"]
r2_read_structure = config["r2_read_structure"]

units_wd = '/'.join([runid, config["units"]])

# 
#validate(config, schema="../../config/schemas/config.schema.yaml")
if Path(units_wd).is_file():
  units = (pd.read_csv(units_wd, sep="\t", dtype={"sample_name": str, "unit_name": str})
      .set_index(["sample_name", "unit_name"], drop=False)
      .sort_index()
  )
  validate(units, schema="../../config/schemas/units.schema.yaml")
  samples_dict = dict(zip(units['sample_name'], zip(units['fq'])))  #, units['fq2'])))
  idkeys = list(samples_dict.keys())
else: 
  #idkeys = list()
  #units= str()
  lsda = os.listdir('/'.join([runid, "results/bcl2fq/"]))
  #lsids = [word for word in lsda if word.endswith("L001_R1_001.fastq.gz") and 'Undetermined' not in word]
  #idkeys = ['_'.join(x.split('_')[0:2]) for x in lsids]
  # create units.tsv#
  lsidsL12 = [word for word in lsda if word.endswith("R1_001.fastq.gz") and 'Undetermined' not in word]
  idkeysL12 = ['_'.join(x.split('_')[0:2]) for x in lsidsL12]
  lanes = [re.findall('L00[12]', x)[0] for x in lsidsL12]
  fq = ['/'.join([runid, 'results/bcl2fq', fqgz]) for fqgz in lsidsL12]
  df = pd.DataFrame({'sample_name': idkeysL12, 'unit_name': lanes, 'fq': fq})
  df.to_csv('/'.join([runid, 'units.tsv']), sep= '\t')
  units = df
  validate(units, schema="../../config/schemas/units.schema.yaml")
  samples_dict = dict(zip(units['sample_name'], zip(units['fq'])))  #, units['fq2'])))
  idkeys = list(samples_dict.keys())
  

# define working directory and samples for rule all:
#pth = Path(units['fq'].iloc[0])
#fqDir = str(pth.parent)

pr = Path('/'.join([runid, 'results/']))
pr.mkdir(parents= True, exist_ok= True)
pl = Path('/'.join([runid, 'logs/']))
pl.mkdir(parents= True, exist_ok= True)


# set wildcard constraints on {sample}
common_constraint = "[0-9A-Za-z\-\_]+[^\/L][^\/umi][^\/chr7][^\/SE]"  #[^/]"

#[^\/mrkdup]

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

onstart:
    shell("gitrev=$(git rev-parse HEAD) && echo  \"--------------------------------------------- \n Running RNAseq fusion workflow, with \n git branch $gitrev \n log file in \n {log} \n & run data in \n {runid} \n  ---------------------------------------------\" 2>&1 | tee {log}") 

onsuccess:
    shell("cp -v {log} {runid}/logs/")

onerror:
    shell("cp -v {log} {runid}/logs/")


