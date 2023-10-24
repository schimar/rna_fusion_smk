# RNA gene fusion & exon skipping - snakemake workflow 


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------

## Prerequisites before running the workflow

# ----------------------
### clone the workflow:
```
git clone tpgit:/home/admin_bio/projects/rna_fusion_quant.git
```

# ----------------------
### install dependencies:
go to the working directory 
(in my case, it was cloned into the folder smk/ in my home directory)
```
cd ~/smk/rna_fusion_quant/workflow/
git pull
```

and install the dependencies for this workflow. For this, you need to have mamba installed (see the [mamba documentation](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install) and [this handy page](https://www.imranabdullah.com/2021-08-21/Conda-and-Mamba-Commands-for-Managing-Virtual-Environments)).

```
mamba env create -f evs/s7_fgbio2.yaml
```
this will create an environment called "fgbio2"; activate it with the following:
```
mamba activate fgbio2
```

# ----------------------
### egfr-v3-determiner

Further, we need this tool for the EGFR variant III detection (install it in a location of your choosing):
```
git clone https://github.com/yhoogstrate/egfr-v3-determiner.git
cd egfr-v3-determiner
python setup.py install 
# and test it with:
nosetests tests/*.py
```

# ----------------------
### Resources

Finally, make sure you have all of the necessary resource files copied into the resources/ folder 
```
# in ~/smk/rna_fusion_quant/workflow/
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
    - # bcldir (location of NovaSeq run folder)
    - # sashdir (location of SampleSheet.csv, if not in bcldir)


### 2) go to working directory
``` 
cd ~/smk/rna_fusion_quant/workflow/
```

### 3) activate mamba environment
```
mamba activate fgbio2
```

### 4) perform dry-run 
```
smk -npr --config runID= units=../config/units_<tmp>.tsv
# e.g.
smk -npr --config runID=/mnt/sda/rnaSeq/runs/230720 units=../config/units_230720.tsv
```

### 5) run workflow
```
smk -j<nthreads> --config runID=<PATH/TO/runid> units=../config/units_<tmp>.tsv
# e.g.
smk -j50 --config runID=/mnt/sda/rnaSeq/runs/230720 units=../config/units_230720.tsv
```


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


## other useful features

### get the directed acyclic graph (DAG) as pdf:
```
smk -npr --config runID=<PATH/TO/runid> units=../config/units_<tmp>.tsv --rulegraph | dot -Tpdf > dag.pdf

# e.g. 
smk -npr --config runID=/mnt/sda/rnaSeq/runs/230720 units=../config/units_230720.tsv --rulegraph | dot -Tpdf > dag.pdf

```


### if you had to cancel a run or if there was an error, append ``--rerun-incomplete`` to the respective smk command 






