# RNA gene fusion & exon skipping - snakemake workflow 


# -----------------------------------------------------------------

## install dependencies:
go to the working directory with 
```
cd ~/smk/rna_fusion_quant/workflow/
```

and install the dependencies for this workflow. For this, you need to have mamba installed (see the [mamba documentation](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install) and [this handy page](https://www.imranabdullah.com/2021-08-21/Conda-and-Mamba-Commands-for-Managing-Virtual-Environments)).

```
mamba env create -f evs/s7_fgbio2.yaml
```
this will create an environment called "fgbio2", activate it with the following:
```
mamba activate fgbio2
```

Further, we need this tool for the EGFR variant III detection (install it in a location of your choosing):
```
git clone https://github.com/yhoogstrate/egfr-v3-determiner.git
cd egfr-v3-determiner
python setup.py install 
# and test it with:
nosetests tests/*.py
```

## prepare the ``units<tmp>.tsv`` file: 

# sashdir
# bcldir
# rundir (-> runid) 
## in the working directory:
```
# e.g. 
~/smk/rna_fusion_quant/workflow/
```

### get the directed acyclic graph (DAG) as pdf:
```
smk -npr --config runID=<PATH/TO/runid> units=../config/units_<tmp>.tsv --rulegraph | dot -Tpdf > dag.pdf

# e.g. 
smk -npr --config runID=/mnt/sda/rnaSeq/runs/230720 units=../config/units_230720.tsv --rulegraph | dot -Tpdf > dag.pdf

```


dry-run: 
```
smk -npr --config runID=<PATH/TO/runid> units=../config/units_<tmp>.tsv
# e.g.
smk -npr --config runID=/mnt/sda/rnaSeq/runs/230720 units=../config/units_230720.tsv
```

run: 
```
smk -j<nthreads> --config runID=<PATH/TO/runid> units=../config/units_<tmp>.tsv
# e.g.
smk -j50 --config runID=/mnt/sda/rnaSeq/runs/230720 units=../config/units_230720.tsv
```

### addtl notes:

- if you had to cancel a run, append ``--rerun-incomplete`` to the respective smk command 






