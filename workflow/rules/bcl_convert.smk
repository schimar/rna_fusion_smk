# -----------------------------------------------------------------------------
# Base calling: bcl-convert (single pass, --no-lane-splitting)
#
# - bcl-convert demuxes ALL samples in the SampleSheet (incl. non-WTS), but
#   this workflow declares outputs ONLY for the WTS samples (wts_samples).
# - After conversion, ALL produced fastqs are renamed/moved out of the fresh
#   tmp dir into the final bcl2fq/ directory (so non-WTS demux is preserved
#   on disk for the future full workflow even though it is not tracked here).
# - Final naming: {sample}.R{1,2}_001.fastq.gz
#   (bcl-convert emits {sample}_S<index>_R{1,2}_001.fastq.gz; the trailing
#   _S<index>_R<read> is converted to .R<read>, keeping {sample} intact —
#   including any _S<n> that is part of the Sample_ID).
# -----------------------------------------------------------------------------

rule bcl_convert:
    input:
        sheet= bcldir + "/SampleSheet_rna.csv",
    output:
        expand("{runid}/results/bcl2fq/{sample}.{read}_001.fastq.gz",
               runid= runid, sample= wts_samples, read=["R1", "R2"]),
    params:
        bcldir= bcldir,
        outdir= runid + "/results/bcl2fq/",
        bcl_tmp= runid + "/results/bcl2fq/.tmp_bcl",
    log:
        runid + "/logs/bcl_convert.log",
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
            --no-lane-splitting true \
            --bcl-sampleproject-subdirectories false \
            --sample-name-column-enabled false \
            --bcl-num-parallel-tiles 32 \
            --bcl-num-conversion-threads {threads} \
            > {log} 2>&1

        # Rename ALL demuxed fastqs into final dir:
        #   SAMPLE_S<index>_R<read>_001.fastq.gz  ->  SAMPLE.R<read>_001.fastq.gz
        python3 scripts/rename_bcl_fastqs.py {params.bcl_tmp} {params.outdir}

        # remove the tmp dir (Reports etc. are regenerated on demand)
        rm -rf {params.bcl_tmp}
        """


# -----------------------------------------------------------------------------
# Read merging for QC (bbmerge)
# Consumes the demuxed gz directly (no cat_lanes, no lane splitting).
# -----------------------------------------------------------------------------
rule bbmerge_fqs:
    input:
        fq1 = "{runid}/results/bcl2fq/{sample}.R1_001.fastq.gz",
        fq2 = "{runid}/results/bcl2fq/{sample}.R2_001.fastq.gz",
    output:
        out   = temp("{runid}/results/bbmerge/{sample}.fq.gz"),
        outu1 = temp("{runid}/results/bbmerge/outu/{sample}.1.fq.gz"),
        outu2 = temp("{runid}/results/bbmerge/outu/{sample}.2.fq.gz"),
        hist  = "{runid}/results/bbmerge/{sample}.hist.txt",
    wildcard_constraints:
        sample = common_constraint,
    log:
        "{runid}/logs/bbmerge/{sample}.log",
    priority: 1
    threads: 12
    shell:
        """
        bbmerge-auto.sh -Xmx60g in1={input.fq1} in2={input.fq2} \
            outm={output.out} outu1={output.outu1} outu2={output.outu2} \
            ihist={output.hist} ecct extend2=20 iterations=5 > {log} 2>&1
        """