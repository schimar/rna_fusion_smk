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

# ----------------------
### install dependencies:
go to the working directory 
(in my case, it was cloned into the folder smk/ in my home directory)
```
cd ~/smk/rna_fusion/workflow/
git pull
```

and install the dependencies for this workflow. For this, you need to have mamba installed (see the [mamba documentation](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install) and [this handy page](https://www.imranabdullah.com/2021-08-21/Conda-and-Mamba-Commands-for-Managing-Virtual-Environments)).

```
mamba env create -f envs/envs.yaml
```
this will create an environment called "fgbio2"; activate it with the following:
```
mamba activate env_rna
```


# ----------------------
### Resources

Finally, make sure you have all of the necessary resource files copied into the resources/ folder 
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

NOTE: currently, we still have to run the bcl2fq rule separately, before running the rest of the workflow (I've made some progress on this, but not quite there yet). 
For this, you need to comment out everything after the bcl2fq file in rule all (in the main snakefile, everything after line 45). Once bcl2fq is done, uncomment the same ones again and you're golden for the second run!


### 1) defne the paths:

(i.e. where is your raw data (the *.bcl files) and where do you want to write the data to?) make sure you have the following info:

    - runid (where to write to - consider writing to local ssd) 
    - bcldir (location of NovaSeq run folder)
    - SampleSheet_rna.csv in bcldir


### 2) go to working directory
``` 
cd ~/smk/rna_fusion/workflow/
```

### 3) activate mamba environment
```
mamba activate env_rna
```

### 4) perform dry-run 
```
smk -npr --config runid=<output_path> bcldir=/mnt/{illumina,routine}/<run_folder>/
# e.g.
smk -npr --config runid=/mnt/sda/rnaSeq/runs/231025_rerun bcldir=/mnt/illumina/231025_A01272_0063_AHK5CGDRX3/
```

### 5) run workflow
```
smk -j<nthreads> --config runID=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/>
# e.g.
smk -j50 --config runID=/mnt/sda/rnaSeq/runs/231025 bcldir=/mnt/illumina/231025_A01272_0063_AHK5CGDRX3/
```


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


## other useful features

### get the directed acyclic graph (DAG) as pdf:
```
smk -npr --config runID=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> --rulegraph | dot -Tpdf > dag.pdf

# e.g. 
smk -npr --config runID=/mnt/sda/rnaSeq/runs/230720 bcldir=/mnt/illumina/230720_A01272_0063_AHK5CGDRX3/ --rulegraph | dot -Tpdf > dag.pdf

```


### if you had to cancel a run or if there was an error, append ``--rerun-incomplete`` to the respective smk command 






