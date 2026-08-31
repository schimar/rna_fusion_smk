# -----------------------------------------------------------------------------
# Read mapping with STAR (whole-transcriptome / fusion-ready alignment).
#
# Consumes the demuxed gzipped fastqs directly (bcl_convert output:
# {runid}/results/fastq/{sample}.R{1,2}_001.fastq.gz). Reads are gzipped, so
# STAR is given --readFilesCommand zcat.
# -----------------------------------------------------------------------------

rule map_star:
    input:
        fq1 = "{runid}/results/fastq/{sample}.R1_001.fastq.gz",
        fq2 = "{runid}/results/fastq/{sample}.R2_001.fastq.gz",
        idx = "resources/star_genome",
        gtf = "resources/genome.gtf",
    output:
        aln = temp("{runid}/results/bam/{sample}.star.bam"),
        sj  = "{runid}/results/bam/{sample}.star/SJ.out.tab",
    params:
        smpl= "{sample}",
        lib = "Library1",
        pu  = "Unit1",
        pl  = "Illumina",
        extra=("--readFilesCommand zcat "
               "--chimSegmentMin 12 --chimOutType WithinBAM "
               "--readFilesSAMattrKeep All --quantMode GeneCounts "
               "--limitBAMsortRAM 110000000000 "
               "--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 "
               "--chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 "
               "--alignSJDBoverhangMin 10 --alignMatesGapMax 100000 "
               "--alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 "
               "--chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 "
               "--chimMultimapNmax 20 --chimNonchimScoreDropMin 10 "
               "--peOverlapNbasesMin 12 --peOverlapMMp 0.1 "
               "--alignInsertionFlush Right "
               "--alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30"),
    log:
        "{runid}/results/bam/{sample}.map_star.log",
    threads: 24
    wildcard_constraints:
        sample = common_constraint,
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} {params.extra} \
            --outFileNamePrefix {runid}/results/bam/{wildcards.sample}.star/ \
            --outStd BAM_SortedByCoordinate \
            --outSAMattrRGline ID:{rgid} SM:{params.smpl} LB:{params.pl} PU:{params.pu} \
            > {output.aln} 2> {log}
        """


rule star_index_dup:
    input:
        "{runid}/results/bam/{sample}.star.bam",
    output:
        "{runid}/results/bam/{sample}.star.bam.bai",
    wildcard_constraints:
        sample = common_constraint,
    log:
        "{runid}/results/bam/{sample}.star_index_dup.log",
    shell:
        """
        samtools index {input} > {log} 2>&1
        """


rule star_markdup:
    input:
        "{runid}/results/bam/{sample}.star.bam",
    output:
        "{runid}/results/bam/{sample}.bam",
    wildcard_constraints:
        sample = common_constraint,
    priority: 20
    threads: 8
    log:
        "{runid}/results/bam/{sample}.star_markdup.log",
    shell:
        """
        sambamba markdup -t {threads} --overflow-list-size 2000000 \
            {input} {output} > {log} 2>&1
        """