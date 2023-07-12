# Bacterial shotgun seq - assembly & annotation - snakemake workflow 


# -----------------------------------------------------------------

dry-run:  
```
smk -npr --config runID=20230125_RNASeq fqDir="/mnt/illumina/230209_A01272_0035_BHTVFGDRX2/Fastq

```

run:  
```
smk -j48 --config runID=20230125_RNASeq fqDir="/mnt/illumina/230209_A01272_0035_BHTVFGDRX2/Fastq
```



