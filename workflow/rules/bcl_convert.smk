# -----------------------------------------------------------------------------
# Base calling: bcl-convert (single pass, --no-lane-splitting)
#
# - bcl-convert demuxes ALL samples in the SampleSheet (incl. non-WTS), but
#   this workflow declares outputs ONLY for the WTS samples (wts_samples).
# - After conversion, ALL produced fastqs are renamed/moved out of the fresh
#   tmp dir into the final fastq/ directory (so non-WTS demux is preserved
#   on disk for the future full workflow even though it is not tracked here).
# - Final naming: {sample}.R{1,2}_001.fastq.gz
#   (bcl-convert emits {sample}_S<index>_R{1,2}_001.fastq.gz; the trailing
#   _S<index>_R<read> is converted to .R<read>, keeping {sample} intact —
#   including any _S<n> that is part of the Sample_ID).
# -----------------------------------------------------------------------------

rule bcl_convert:
    input:
        sheet=sample_sheet,
        samples=samples_tsv,
    output:
        expand("{runid}/results/fastq/{sample}.{read}_001.fastq.gz",
               runid= runid, sample= wts_samples, read=["R1", "R2"]),
    params:
        bcldir= bcldir,
        outdir= runid + "/results/fastq/",
        bcl_tmp= runid + "/results/fastq/.tmp_bcl",
        no_lane_splitting=no_lane_splitting,
    log:
        runid + "/results/fastq/bcl_convert.log",
    threads: 56
    shell:
        """
        set -eo pipefail
        mkdir -p "$(dirname {log})"

        # bcl-convert 4.x requires the output directory to NOT already exist
        rm -rf {params.bcl_tmp}
        mkdir -p "$(dirname {params.bcl_tmp})"

        bcl-convert \
            --bcl-input-directory {params.bcldir} \
            --sample-sheet {input.sheet} \
            --output-directory {params.bcl_tmp} \
            --no-lane-splitting {params.no_lane_splitting} \
            --bcl-sampleproject-subdirectories false \
            --sample-name-column-enabled false \
            --bcl-num-parallel-tiles 32 \
            --bcl-num-conversion-threads {threads} \
            > {log} 2>&1

        # Rename ALL demuxed fastqs into final dir:
        #   SAMPLE_S<index>[_L<lane>]_R<read>_001.fastq.gz
        #       -> SAMPLE[_L<lane>].R<read>_001.fastq.gz
        python3 scripts/rename_bcl_fastqs.py {params.bcl_tmp} {params.outdir}

        # keep demux metadata alongside the fastqs (checkpoint convention)
        cp -r {params.bcl_tmp}/Reports {params.bcl_tmp}/Stats {params.outdir} 2>/dev/null || true
        cp {params.bcldir}/RunInfo.xml {params.bcldir}/RunParameters.xml {params.outdir} 2>/dev/null || true

        # remove the tmp dir
        rm -rf {params.bcl_tmp}
        """


# -----------------------------------------------------------------------------
# Read merging for QC (bbmerge)
# Consumes the demuxed gz directly (no cat_lanes, no lane splitting).
# -----------------------------------------------------------------------------
rule bbmerge_fqs:
    input:
        fq1 = "{runid}/results/fastq/{sample}.R1_001.fastq.gz",
        fq2 = "{runid}/results/fastq/{sample}.R2_001.fastq.gz",
    output:
        out   = temp("{runid}/results/quality_control/fastq/{sample}.bbmerge.fq.gz"),
        outu1 = temp("{runid}/results/quality_control/fastq/{sample}.bbmerge.unmerged.1.fq.gz"),
        outu2 = temp("{runid}/results/quality_control/fastq/{sample}.bbmerge.unmerged.2.fq.gz"),
        hist  = "{runid}/results/quality_control/fastq/{sample}.bbmerge.hist.txt",
    wildcard_constraints:
        sample = common_constraint,
    log:
        "{runid}/results/quality_control/fastq/{sample}.bbmerge_fqs.log",
    priority: 1
    threads: 12
    shell:
        """
        bbmerge-auto.sh -Xmx60g in1={input.fq1} in2={input.fq2} \
            outm={output.out} outu1={output.outu1} outu2={output.outu2} \
            ihist={output.hist} ecct extend2=20 iterations=5 > {log} 2>&1
        """