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
  `SampleSheet.csv`, or a separate SampleSheet path supplied as `sample_sheet`
- `runid` (optional): workflow output directory; defaults to the `bcldir`
  basename (consider an explicit local-SSD path for larger runs)
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

### Reference genome

The reference genome is selected in `config/config.yaml` via the `ref` key
(file stem). The default is
`GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18`, so
`resources/` must contain, named after that stem:

- `<ref>.fasta` — reference sequence
- `<ref>.gtf` — annotation
- `<ref>.txt` — chromosome names/sizes (used by `bedtools sort -faidx` / `-g`)
- `star_<ref>/` — STAR index (built by the workflow's `star_index` rule)

Arriba selects its bundled, version-matched blacklist from `genome_build`.
To run a different genome, provide the files under the new stem and set `ref:`.

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
- `final_dest` (optional): destination directory for delivered results; defaults to
  `/home/schilling_m1/smb/Analyses/00_Tests/vc_rna/<bcldir basename>`
- `rgid` (optional): read-group ID; defaults to the `bcldir` basename
- `sample_sheet` (optional): path to the SampleSheet; defaults to
  `<bcldir>/SampleSheet.csv`
- `lane_splitting` (optional): `false` by default; set `true` only for a
  lane-split bcl-convert run

### 2) perform a dry-run

```
cd <repo>/workflow
../sub-smk-job.sh -n --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> sample_sheet=<PATH/TO/SampleSheet.csv> final_dest=<PATH/TO/final_dest>
```

### 3) run the workflow

```
../sub-smk-job.sh -j50 --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> sample_sheet=<PATH/TO/SampleSheet.csv> final_dest=<PATH/TO/final_dest>
```

The workflow finishes when `{final_dest}/run.done` exists. Results are rsync'd
to `final_dest`.

## Sample definition

- The `sample_sheet` config value, or `<bcldir>/SampleSheet.csv` by default,
  is the source of truth for the sample list (both v1 and v2 sheet layouts are
  supported).
- A sample is treated as a WTS sample iff its full `Sample_ID` contains `WTS`.
  The rest of the ID (patient identifiers, `_S<n>` suffix, ...) is not
  interpreted.
- `bcl-convert` still demultiplexes all samples of the sheet; only the WTS
  fastqs become workflow-managed outputs, so pooled non-WTS fastqs remain on
  disk for other (future) workflows.
- A derived `samples.tsv` (sample_id + is_wts) is written into `runid/` for
  transparency.

## Results layout (checkpoint files)

Results are written in a **checkpoint-file layout** (aligned with the varCAD
workflow): one top-level directory per file type, the sample identity lives in
the **filename**, QC mirrors the type of the file it was derived from, and each
rule's log is written **next to its output** as `{output}.{rule}.log`
(universal rules like `bcl_convert`/`multiqc` get `{dir}/{rule}.log`;
reference-build logs stay colocated in `resources/`).

```
{runid}/results/
  fastq/     {sample}.R{1,2}_001.fastq.gz (+ Reports/ Stats/ InterOp/ RunInfo.xml RunParameters.xml)
  bam/       {sample}.bam | {sample}.bam.bai | {sample}.star/SJ.out.tab
  fusions/   {sample}.fusions.tsv | {sample}.fusions.discarded.tsv (+ other arriba artifacts)
  splicing/  {sample}.SE.MATS.JC.txt | {sample}.egfr_v3.out
  quality_control/
    fastq/   {sample}.bbmerge.hist.txt | {sample}.clump.{opt,pcr}.stats.txt | multiqc_report.html
    bam/     {sample}.bcov.tsv | {sample}.n10.tsv
```

The alignment checkpoint is a **BAM** for now (STAR's native output; a CRAM
checkpoint may be added later). The final `rsync` simply mirrors
`{runid}/results/` into `{final_dest}/`, so this layout is
preserved in the delivered results.

## other useful features

### get the rulegraph (preferred over the full DAG) as png:
```
../sub-smk-job.sh -n --config runid=<PATH/TO/runid> bcldir=<PATH/TO/bcldir/> analysis_path=analysis_rna --rulegraph --quiet | dot -Tpng > rulegraph.png
```

### if you had to cancel a run or if there was an error, append ``--rerun-incomplete`` to the respective command






