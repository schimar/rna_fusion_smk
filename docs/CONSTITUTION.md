# WTS-Init — Design & Constitutional Specification (SPEC)

Status: **ACCEPTED** — living document; changes are additive and reviewed.
Branch target: `wts-init` (off clean `master` at `1538fb8`).
Review/approval: this SPEC must be reviewed and approved before implementation begins; each implementation commit must reference the SPEC requirement(s) it satisfies.

---

## 0. Purpose

Refactor the `rna_fusion` Snakemake RNA workflow into a **whole-transcriptome-sequencing (WTS) init** pipeline that runs entirely inside a single central Docker image (`ukwgenommedizin/rna_fusion:1.0.1`), drops the fgbio/UMI path, and emits raw fusion / alternative-splicing outputs without the post-hoc filtering rules.

The SPEC exists to give a precise, verifiable contract for the change — the antithesis of ad-hoc "hand-wavey" edits that cause endless debugging rounds.

## 1. Scope

### 1.1 In scope
- Instantiate the workflow under a **single central Docker image** containing Snakemake 9.19.0 + all required tools.
- Provide a launcher script (`sub-smk-job.sh`) that runs Snakemake **inside** that image.
- Read demultiplexed FASTQ via **`bcl-convert --no-lane-splitting true`** (drop `bcl2fastq`), with an in-rule **rename** to `{sample}.R{1,2}_001.fastq.gz`.
- Drop the fgbio/UMI path; move the STAR rules into `mapping.smk` and `map_star` consumes the demuxed **gzipped** R1/R2 directly (`--readFilesCommand zcat`).
- Emit **raw** Arriba, rMATS, and EGFRv3 outputs; keep the post-hoc filter rules defined but **inert** (not in the DAG).
- Keep QC outputs in `run.done`.
- Ensure compatibility with **Snakemake 9.19.0** (drop the top-level `report:` directive; keep `min_version("7.0.0")`).

### 1.2 Explicitly out of scope (non-goals — do not re-litigate)
- No per-rule `docker run`; no Docker-outside-Docker / orchestrator-calls-docker.
- No `--use-conda` / `conda:` directives at runtime; no per-tool image split.
- No Snakemake executor plugins; no SLURM/cluster profile. Local executor only.
- **No wrapper → shell rewrites.** Wrapper rules invoke tools by name; those tools exist in the image and run as-is.
- No deletion of dead rules/scripts in this branch (cleanup is a follow-up branch).
- No change to `master`; work isolated to `wts-init`.

## 2. Environment & Image Contract

- The single image `ukwgenommedizin/rna_fusion:1.0.1` (fallback to `1.0.0`) provides, on `PATH`:
  - snakemake 9.19.0, python3 3.12.14
  - STAR 2.7.11b, arriba 2.5.1, bbmap 39.01 (bbmerge.sh/clumpify.sh)
  - fastqc, bedtools, sambamba, bwa, bwa-mem2, samtools
  - multiqc 1.35, rmats.py, egfr-v3-determiner 0.7.4, bcl-convert 4.5.4, fastp
  - `fgbio.jar` and `picard.jar` under `/opt/tools/jars`
- Absent from image (so must not be required): `bcl2fastq`, `cutadapt`, `picard`-on-PATH, `fgbio`-on-PATH (jars only).
- Image default is overridable via env var `RNA_FUSION_IMAGE`.

## 3. Execution Model

### 3.1 Launcher — `sub-smk-job.sh` (repo root, executable)
A single `docker run` executes Snakemake inside the central image and passes CLI args (including `--config`) straight through.

Acceptance: `./sub-smk-job.sh -n --config runid=... bcldir=... final_dest=...` performs a dry-run; no `--use-conda`, no per-rule docker, no SLURM.

### 3.2 Mounts & user
- Mount `/home`, `/mnt`, and the repo path (read/write for run outputs).
- Run as invoking user: `-u $(id -u):$(id -g)` so outputs are owned by the caller.

## 4. Demultiplexing & read handling

### 4.1 `bcl_convert` rule (rewrite; keep; `bcl2fq` dropped)
- Use `bcl-convert 4.5.4` with:
  - `--bcl-input-directory {bcldir}`
  - `--sample-sheet {sample_sheet}` (defaults to `{bcldir}/SampleSheet.csv`)
  - `--output-directory <fresh tmp dir>` (must not pre-exist; bcl-convert 4.x requirement)
  - `--no-lane-splitting true`
  - `--bcl-sampleproject-subdirectories false`
  - `--bcl-num-parallel-tiles ...` (tune)
- After bcl-convert returns, **rename loop** (same-filesystem `mv`, metadata-only): `tmp/{sample}_S*_R{1|2}_001.fastq.gz` → `{runid}/results/bcl2fq/{sample}.R{1|2}_001.fastq.gz`. Exactly one merged file per sample+direction is expected under no-lane-splitting.
- Drop `bcl2fastq` rule, `bcl2fastq.yaml` env, and `cat_lanes`.

Acceptance: demuxed files exist at `{runid}/results/bcl2fq/{sample}.R1_001.fastq.gz` and `.R2_001.fastq.gz`.

### 4.2 Sample definition & selection (replaces `units.tsv`)

- **Source of truth = SampleSheet** (`{sample_sheet}`, defaulting to
  `{bcldir}/SampleSheet.csv`), read at **parse time**.
- Support **both v1 and v2** sheet layouts:
  - v1: header row contains `Sample_ID` at column 0.
  - v2 (NovaSeqX+): locate `Sample_ID` by header name (may be any column).
  - Robust approach (same as vrc reference): scan lines for the header row containing `Sample_ID`, then collect non-empty values in that column.
- **Routing key = full `Sample_ID` substring match.** Sample is **WTS** iff `"WTS" in sample_id`. Do **not** interpret any other part of the ID (patient, project, `_S<n>`, etc.) — those are opaque.
- **`Sample_ID` is used whole as the single identifier everywhere** (wildcard `{sample}`, fastq paths, downstream outputs). Example: `B-WTS26-0031_Pat26-3689_26-004199-B-I-I_S1` (kept as one ID; `_S<n>` is the internal-control suffix and is retained).
- `wts_samples` = the list of matching IDs; `idkeys` = `wts_samples`; all downstream `expand(...)` and wildcard rules use it.
- The derived sample table is written to `<runid>/samples.tsv` for transparency/validation (derived artifact, not source of truth).

Acceptance: `snakemake -n` builds the full DAG (bcl-convert → mapping → QC → fusions/altsplice) in **one** invocation, targeting only `wts_samples`, with no two-phase bcl-convert hack.

### 4.3 `bcl_convert` — WTS-only DAG, all-sample physical demux
- `rule bcl_convert` declares **outputs only for `wts_samples`** (predicted paths):
  ```
  output: expand("{runid}/results/bcl2fq/{sample}.R{read}_001.fastq.gz", sample=wts_samples, read=[1,2])
  ```
- **Physically**, `bcl-convert` still demuxes **all** samples in the sheet (`--no-lane-splitting true`), and the rename loop moves **all** produced fastqs out of the fresh tmp dir to `{runid}/results/bcl2fq/` — so non-WTS demux is preserved on disk (consumed later by the large workflow), even though Snakemake does not track it.
- Declared outputs drive job creation/ordering; the DAG is fully known at parse time; downstream rules reference the predicted fastq paths directly → **single-run execution**.

### 4.3.1 Lane-split test runs
- `lane_splitting: false` remains the default and yields one workflow unit per
  full `Sample_ID`.
- With `lane_splitting: true`, each `Lane` value becomes a separate workflow
  unit named `<full Sample_ID>_L###`. The full opaque `Sample_ID`, including
  any existing `_L...` portion, is preserved; the appended suffix identifies
  the bcl-convert sequencer lane. Lanes are not merged.

### 4.4 `common.smk` units parsing
- Replace `detect_lanes_and_create_units` (file-scan → parse-time table): remove `units.tsv` dependency.
- Remove `units` key from `config.yaml`; `units.schema.yaml` stays on disk but unused.
- Warn if a `Sample_ID` contains a `.` (would break later path parsing) — README note + runtime warning.

## 5. Dropping UMI; adding `mapping.smk`

- New `workflow/rules/mapping.smk` contains the 3 STAR rules moved from `umi.smk`: `map_star`, `star_index_dup`, `star_markdup`.
- `map_star` inputs become the demuxed gz: `{runid}/results/bcl2fq/{sample}.R{1,2}_001.fastq.gz`, and STAR `extra` gains `--readFilesCommand zcat`.
- `snakefile`: remove `include: "rules/umi.smk"`, add `include: "rules/mapping.smk"`. `umi.smk` stays on disk (un-included).
- fgbio-specific QC targets trimmed from active DAG (`fastqc_cons`, `multiqc_cons`).

Acceptance: `snakemake -n` produces no fgbio/UMI/consensus rules; STAR runs read the gz via zcat.

## 6. Outputs — raw deliverables, filters inert

### 6.1 Arriba
- Active deliverable: raw `fusions.tsv` (+ `fusions.discarded.tsv`).
- Rules `get_clinFuse`, `filter_arriba`, `getExons_arriba`, `gtf2exon_h5` stay **defined but unreachable** (not in `rule all`).

### 6.2 rMATS
- Active deliverable: raw per-chrom `{sample}/{chrom}/SE.MATS.JC.txt`.
- Rules `clinEx`, `exonfltr`, `cat_exonSkippers`, `exonFinalOut` stay **defined but unreachable**.

### 6.3 EGFR
- `{sample}.egfr_v3.out` unchanged (already a leaf output).

### 6.4 `rule all` / `out.smk` (rsync → `run.done`)
- Repoint to raw deliverables + kept QC: `fastqc/multiqc`, `bcov`, `n10cov`.
- Drop `.fltrd.ex` / `.exonSkip` references.

## 7. Snakemake 9.19 compatibility
- Drop the top-level `report: "report/workflow.rst"` directive (removed in SMK8+).
- Keep `min_version("7.0.0")`.
- Ignore `config/slurm/*` (not used; local executor).

## 8. Validation / Definition of Done
1. `git switch -c wts-init` (done) — first commit is THIS SPEC.
2. `./sub-smk-job.sh -n --config ...` dry-run succeeds with no `--use-conda` warnings and no fgbio/UMI rules in the DAG.
3. A small (1–2 sample) end-to-end run produces: demuxed `{sample}.R{1,2}_001.fastq.gz`, STAR BAM + `SJ.out.tab` + `mrkdup.bam`, arriba `fusions.tsv`(+`discarded.tsv`), rMATS `SE.MATS.JC.txt`, `egfr_v3.out`, QC outputs, and `{final_dest}/run.done`.
4. No `bcl2fastq`, `cat_lanes`, fgbio, or filter rules in the executed DAG.
5. Incremental commits, each referencing its SPEC section(s).

## 9. Open Decision Register (resolved where marked)
- D1 `RNA_FUSION_IMAGE` default → `1.0.1` (RESOLVED).
- D2 bcl-convert parallel-tiles value → tune, default from reference (~32) (PENDING measure at runtime).
- D3 Dot-in-sample-name handling → warn + document (RESOLVED).
- D4 QC set in `run.done` → fastqc/multiqc, bcov, n10cov (RESOLVED).

## 10. Amendment A — checkpoint-file results layout (SUPERSEDES §6.1–6.4 paths)

All results move from tool-centric dirs to a **checkpoint-file layout**
(aligned with the varCAD workflow): top-level dir = checkpoint file type,
sample identity lives in the filename, QC mirrors its source file type under
`quality_control/{type}/`, and logs are **colocated with their output** as
`{output_file}.{rule}.log` (universal rules: `{dir}/{rule}.log`; reference
builds: colocated in `resources/`).

| Checkpoint | Path |
|---|---|
| fastq | `results/fastq/{sample}.R{1,2}_001.fastq.gz` (+ `Reports/`, `Stats/`, `InterOp/`, `RunInfo.xml`, `RunParameters.xml`) |
| alignment | `results/bam/{sample}.bam` + `.bam.bai` (STAR aux in `bam/{sample}.star/`) |
| fusions | `results/fusions/{sample}.fusions.tsv` / `.discarded.tsv` (+ other arriba artifacts) |
| splicing | `results/splicing/{sample}.SE.MATS.JC.txt` / `{sample}.egfr_v3.out` |
| QC | `results/quality_control/fastq/{sample}.bbmerge.hist.txt`, `{sample}.clump.{opt,pcr}.stats.txt`, `multiqc_report.html`; `results/quality_control/bam/{sample}.bcov.tsv`, `{sample}.n10.tsv` |

- D5 Alignment checkpoint format → **BAM now, CRAM later** (RESOLVED). The raw
  STAR alignment is a temp `bam/{sample}.star.bam`; `bam/{sample}.bam` is the
  markdup checkpoint consumed downstream.
- D6 Arriba `fusions.vcf.gz` location → **`fusions/`** together with the TSVs
  (RESOLVED; arriba artifacts stay together).
- D7 rMATS/EGFR outputs → both under **`splicing/`** (RESOLVED). `rMats` writes
  into a per-sample scratch dir `.rmats_{sample}/` and the result is moved to
  the sample-tagged `splicing/{sample}.SE.MATS.JC.txt` (rmats.py writes a flat
  filename); the placeholder-touch fallback was removed.
- `out.smk` rsync now mirrors `{runid}/results/` → `{final_dest}/`
  (varCAD ExportFolder equivalent), preserving the checkpoint layout.
- Integrity sidecars (`.size`, `.sha256`) — DEFERRED to a later branch.
- Old-path rules (`umi.smk`, `sub_chr` chain, filter rules, rseqc/deeptools)
  are inert and keep their old paths (cleanup branch).

## 11. Amendment B — configurable reference genome

- The reference genome is selected in `config.yaml` via a single `ref` key =
  file stem under `workflow/resources/`. Derived names:
  `{ref}.fasta`, `{ref}.gtf`, `{ref}.txt` (chrom sizes for bedtools),
  `{ref}.dict`, `{ref}.fasta.fai`, and STAR index dir `resources/star_{ref}/`.
- Default (`D8`, RESOLVED):
  `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18`.
- Arriba selects its bundled, version-matched blacklist from the configured
  `genome_build`; no image-internal package path is configured.
- Globals defined once in `common.smk` (`ref`, `genome_fa`, `genome_gtf`,
  `genome_txt`, `genome_dict`, `genome_fai`, `star_idx`); all active rules (`star_index`, `sort_bed`,
  `map_star`, `arriba`, `rMats`, `targetcov_bed`, `rule all`) consume the
  globals. Inert legacy rules keep literal old paths (cleanup branch).
- The old `ref: {species, release, build}` download settings are preserved in
  config as `legacy_ref:` (unused).
- DAG/rulegraph image artifacts (`workflow/dag*.png`, `workflow/rulegraph*.png`)
  are git-ignored; going forward only the **rulegraph** is generated (`D9`,
  RESOLVED — user preference over the full DAG).

## 12. Amendment C — QC (no capture design; transcript-space QC)

### 12.1 Reference pairing (extends §1)
The single `ref` stem in config keeps FASTA/GTF as one pair from
`resources/{ref}.fasta` / `{ref}.gtf`. See §1–2 for the underlying
consistency and no-hardcoding rules (unchanged).

The workflow normalizes the source GTF's contig names against the configured
FASTA before STAR and annotation-derived QC consume it. This supports Ensembl
names such as `1` and `MT` with a `chr1`/`chrM` FASTA.

### 12.2 On-target removal
`sort_bed` / `targetcov_bed` / `cov_n10` are removed from the active DAG
(defined but inert); `config["bed"]` (capture panel) is deprecated/inert.

### 12.3 rRNA BED build guard (extends §3)
rRNA BED generation MUST warn if zero features are found for the in-use GTF.

### 12.4 Globin transcript selection fallback (extends §2, §3)
When building the HBA/HBB subset BED, select transcripts from the in-use
GTF tags in this order: `MANE_Select`/`MANE_Plus_Clinical` first,
`Ensembl_canonical` second, bare `gene_name ∈ {HBA1,HBA2,HBB}` last. The
bare-`gene_name` fallback MUST log loudly — it is the only path that can
select multiple transcripts per gene. All three paths remain tag-derived
from the in-use GTF; none hardcode accessions.

### 12.5 MultiQC content — this round
- In the report: `read_distribution.py` (whole-transcriptome), `read_distribution.py`
  (globin aggregate), rRNA contamination fraction.
- MultiQC natively parses the two RSeQC read_distribution logs. The rRNA
  fraction is a custom overlap check (bedtools coverage / samtools view -c -L),
  not a native RSeQC module — it MUST be surfaced via a MultiQC **custom-content**
  section (JSON/YAML sidecar), not assumed to auto-parse.

### 12.6 Deferred to next increment (per ADR-010)
- BAM-only RSeQC modules: `read_duplication.py`, `read_GC.py`, `bam_stat.py`.
- BED12-dependent trio: `inner_distance.py`, `infer_experiment.py`,
  `junction_annotation.py`, `junction_saturation.py`.
- deeptools `plotCoverage`/`multiBamSummary` coverage-depth (vs derived flat
  exonic BED3).
- This round's multiqc report is intentionally thinner than the eventual
  target; the deferral is a deliberate smaller-commit trade-off, not a gap.

### 12.7 MANE transcript-space read distribution
- Build `resources/mane_{ref}.bed` from every genome-wide transcript tagged
  `MANE_Select` or `MANE_Plus_Clinical`; warn if the in-use GTF has none.
- Run one aggregate `read_distribution.py` report against that BED for each
  sample, alongside the whole-transcriptome and globin reports.
