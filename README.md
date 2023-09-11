# RNA gene fusion & exon skipping - snakemake workflow 


# -----------------------------------------------------------------
## prepare the ``units<tmp>.tsv`` file: 


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



