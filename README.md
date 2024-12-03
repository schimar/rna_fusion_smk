# RNA gene fusion & exon skipping - snakemake workflow 


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------

## Prerequisites before running the workflow

# ----------------------
### clone the workflow:
```
git clone tpgit:/home/admin_bio/projects/rna_fusion.git
```
NOTE: that you need to specify the user and IP adress of the tpgit server if you do not have it in your .ssh/config

# ----------------------
### Resources

Make sure you have all of the necessary resource files copied into the resources/ folder 

```
# in ~/smk/rna_fusion/workflow/
rsync -avzP /mnt/routine/pipelines/shared_resources/rna_v1/* resources/
```


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

## Running the workflow: 

  1) define the paths
  2) go to working directory
  3) activate mamba environment (if not already done)
  4) perform dry-run
  5) run workflow

# -----------------------------------------------------------------------------

NOTE: we have to run the bcl2fq rule separately, before running the rest of the workflow (I've made some progress on this, but not quite there yet). 
For this, you need to comment out everything after the bcl2fq file in rule all (in the main snakefile, everything after line 45). Once bcl2fq is done, uncomment the same ones again and you're golden for the second run!


### 1) defne the paths:

(i.e. where is your raw data (the *.bcl files) and where do you want to write the data to?) make sure you have the following info:

    - runid (where to write to - consider writing to local ssd) 
    - bcldir (location of NovaSeq run folder - the final output will be copied into this folder)
    - SampleSheet_rna.csv in bcldir

NOTE: that if you are running on a panel different than Twist Exome 2.0, you currently need to specify two more parameters, namely <bed> and <analysis_path>. <bed> specifies the target regions (e.g. resources/twist_rna_fusion_target_merged_hg38_annot.bed for the CeGaT panel) and <analysis_path> specifies the name of the output folder, which will be written in <bcldir> (this is only necessary, if you have exome2.0 as well as cegat samples in the same sequencing run). 

### 2) go to working directory
``` 
cd ~/smk/rna_fusion/workflow/
```

### 3) activate mamba environment
```
mamba activate /opt/snakemake
```

### 4) perform dry-run 
```
smk -np --use-conda --conda-prefix /opt/envs/ --conda-frontend mamba --config runid=<output_path> bcldir=/mnt/{illumina,routine}/<run_folder>/
# e.g.
smk -np --use-conda --conda-prefix /opt/envs/ --conda-frontend mamba --config runid=/mnt/sda/rnaSeq/runs/231025_rerun bcldir=/mnt/illumina/231025_A01272_0063_AHK5CGDRX3/
```

### 5) run workflow
```
smk -j<nthreads> --use-conda --conda-prefix /opt/envs/ --conda-frontend mamba --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/>
# e.g.
smk -j50 --use-conda --conda-prefix /opt/envs/ --conda-frontend mamba --config runid=/mnt/sda/rnaSeq/runs/231025 bcldir=/mnt/illumina/231025_A01272_0063_AHK5CGDRX3/
```


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


## other useful features

### get the directed acyclic graph (DAG) as pdf:
```
smk -np --use-conda --conda-prefix /opt/envs/ --conda-frontend mamba --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> --dag --quiet | dot -Tpng > dag.png

# e.g. 
smk -np --use-conda --conda-prefix /opt/envs/ --conda-frontend mamba --config runid=/mnt/sda/rnaSeq/runs/230720 bcldir=/mnt/illumina/230720_A01272_0063_AHK5CGDRX3/ --dag --quiet | dot -Tpng > dag.png

```


### if you had to cancel a run or if there was an error, append ``--rerun-incomplete`` to the respective smk command 






