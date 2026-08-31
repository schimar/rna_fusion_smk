# RNA gene fusion & exon skipping — WTS (whole-transcriptome-sequencing) snakemake workflow


## Running the workflow

Everything runs inside a single central Docker image
(`ukwgenommedizin/rna_fusion:1.0.1`), which contains Snakemake 9.19.0 and all
required bioinformatics tools. There is no conda / mamba activation and no
`--use-conda`; the workflow is launched with the `sub-smk-job.sh` launcher.

The workflow is **single-run**: the SampleSheet is parsed at startup and the
sample set is determined from it (samples whose `Sample_ID` contains `WTS`).
`bcl-convert` demultiplexes the run folder (`--no-lane-splitting`), all WTS
fastqs are renamed to `{sample}.R{1,2}_001.fastq.gz`, and every downstream rule
(mapping, QC, fusion, exon skipping) runs in the same invocation. `bcl2fastq`
is not used.

## Prerequisites before running the workflow

- A NovaSeq run folder (`bcldir`) containing the BCL data and a
  `SampleSheet_rna.csv`
- `runid`: where to write results (consider local SSD)
- The resource/reference files in `workflow/resources/` (genome, indices, bed,
  gene lists). See below.

### clone the workflow:
```
git clone tpgit:/home/admin_bio/projects/rna_fusion.git
```
NOTE: you need to specify the user and IP address of the tpgit server if you do not have it in your `.ssh/config`

### Resources

Make sure you have all of the necessary resource files copied into the
`resources/` folder:

```
# in ~/smk/rna_fusion/workflow/
rsync -avzP /mnt/routine/pipelines/shared_resources/rna_v1/* resources/
```

### Docker image

The image must exist on the host (pull once before running):
```
# (optional) pull the default image
# docker pull ukwgenommedizin/rna_fusion:1.0.1
```
You can override the image via the environment variable `RNA_FUSION_IMAGE`.

### 1) define the paths

You need the following info:

- `runid`: where to write results (consider local SSD)
- `bcldir`: location of the NovaSeq run folder (final output is copied back there)
- `SampleSheet_rna.csv` in `bcldir`

### 2) perform a dry-run

```
cd <repo>/workflow
../sub-smk-job.sh -n --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> analysis_path=analysis_rna
```

### 3) run the workflow

```
../sub-smk-job.sh -j50 --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> analysis_path=analysis_rna
```

The workflow finishes when `{bcldir}/{analysis_path}/run.done` exists
(results are rsync'd back to `bcldir`).

## Sample definition

- The `SampleSheet_rna.csv` in `bcldir` is the source of truth for the sample
  list (both v1 and v2 sheet layouts are supported).
- A sample is treated as a WTS sample iff its full `Sample_ID` contains `WTS`.
  The rest of the ID (patient identifiers, `_S<n>` suffix, ...) is not
  interpreted.
- `bcl-convert` still demultiplexes all samples of the sheet; only the WTS
  fastqs become workflow-managed outputs, so pooled non-WTS fastqs remain on
  disk for other (future) workflows.
- A derived `samples.tsv` (sample_id + is_wts) is written into `runid/` for
  transparency.

## other useful features

### get the directed acyclic graph (DAG) as pdf:
```
../sub-smk-job.sh -n --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> analysis_path=analysis_rna --dag --quiet | dot -Tpng > dag.png
```

### if you had to cancel a run or if there was an error, append ``--rerun-incomplete`` to the respective command






