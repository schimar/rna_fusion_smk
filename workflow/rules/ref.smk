#rule get_genome:
#    output:
#        "resources/genome.fa",
#    log:
#        "logs/get_genome.log",
#    params:
#        species=config["ref"]["species"],
#        datatype="dna",
#        build=config["ref"]["build"],
#        release=config["ref"]["release"],
#    cache: True
#    wrapper:
#        "0.77.0/bio/reference/ensembl-sequence"
#
#
#rule get_annotation:
#    output:
#        "resources/genome.gtf",
#    params:
#        species=config["ref"]["species"],
#        fmt="gtf",
#        build=config["ref"]["build"],
#        release=config["ref"]["release"],
#        flavor="",
#    cache: True
#    log:
#        "logs/get_annotation.log",
#    wrapper:
#        "0.77.0/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        genome_fa,
    output:
        genome_fai,
    log:
        "resources/logs/genome_faidx.log",
    cache: True
    wrapper:
        "0.77.0/bio/samtools/faidx"


rule bwa_index:
    input:
        genome_fa,
    output:
        multiext(genome_fa, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        "resources/logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v1.23.5/bio/bwa-mem2/index"


rule star_index:
    input:
        fasta=genome_fa,
        annotation=genome_gtf,
    output:
        directory(star_idx),
    threads: 8
    params:
      extra=lambda wc, input: f"--sjdbGTFfile {genome_gtf} --sjdbOverhang 100",
    log:
        "resources/logs/star_index.log",
    cache: True
    wrapper:
        "v1.28.0/bio/star/index"  
        #"0.77.0/bio/star/index"



rule create_dict:
    input:
        genome_fa,
    output:
        genome_dict,
    log:
        "resources/logs/create_dict.log",
    params:
        extra="",  # opti.ASM852v1onal: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_gb=8,
    wrapper:
        "v2.1.1/bio/picard/createsequencedictionary"


rule bed_to_interval_list:
    input:
        bed="resources/{bed}_file_UMI_demo_data_hg38.bed",
        dict=genome_dict,
    output:
        "resources/{bed}_file_UMI_demo_data_hg38.interval_list",
    log:
        "resources/logs/{bed}_file_UMI_demo_data_hg38.log",
    params:
        # optional parameters
        extra="--SORT true",  # sort output interval list before writing
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v2.1.1/bio/picard/bedtointervallist"

rule sort_bed:
  input:
    bed = config["bed"],
    chroms = genome_txt
  output:
    "{bed_file_stem}.srtd.bed"
  log: "{bed_file_stem}.sort_bed.log"
  shell:"""
    bedtools sort -faidx {input.chroms} -i {input.bed} > {output} 2> {log}
    """
    #sort -k1,1V -k2,2n -k3,3n {input} > {output}


# -----------------------------------------------------------------------------
# QC reference-build BEDs (SPEC §12; all derived from the in-use {ref}.gtf —
# see CONSTITUTION §12.2/§12.3/§12.4 for the rules). Cacheable.
rule normalize_gtf_contigs:
    input:
        gtf=genome_gtf_source,
        fasta=genome_fa,
    output:
        genome_gtf,
    log:
        f"resources/logs/{ref}.normalized.gtf.normalize_gtf_contigs.log",
    shell:
        "python3 scripts/normalize_gtf_contigs.py -g {input.gtf} -f {input.fasta} -o {output} 2> {log}"


rule gtf2bed_wts:
    input:
        gtf=genome_gtf,
    output:
        "resources/annotation_{ref}.bed",
    log:
        "resources/logs/annotation_{ref}.bed.gtf2bed_wts.log",
    cache: True,
    shell:
        """
        gffread --bed {input.gtf} -o {output} > {log} 2>&1
        """


# HBA/HBB MANE-transcript subset (aggregate BED12; one read_distribution result
# over it, never per-gene — SPEC §12.2-2/§12.3)
rule sub_mane_globin:
    input:
        gtf=genome_gtf,
    output:
        bed="resources/globin_{ref}.bed",
        subset=temp("resources/globin_{ref}.subset.gtf"),
    log:
        "resources/logs/globin_{ref}.bed.sub_mane_globin.log",
    shell:
        """
        python3 scripts/sub_mane_globin.py -g {input.gtf} -o {output.subset} > {log} 2>&1 &&
        gffread --bed {output.subset} -o {output.bed} >> {log} 2>&1
        """


# Genome-wide MANE-tagged transcript subset (aggregate BED12 for clinical
# canonical transcript-space read distribution).
rule mane_bed:
    input:
        gtf=genome_gtf,
    output:
        bed="resources/mane_{ref}.bed",
        subset=temp("resources/mane_{ref}.subset.gtf"),
    log:
        "resources/logs/mane_{ref}.bed.mane_bed.log",
    shell:
        """
        python3 scripts/filter_mane_transcripts.py -g {input.gtf} -o {output.subset} > {log} 2>&1 &&
        gffread --bed {output.subset} -o {output.bed} >> {log} 2>&1
        """


# rRNA loci (flat BED6, overlapping intervals merged; warn-if-zero guard —
# SPEC §12.2-3/§12.3)
rule rRNA_bed:
    input:
        gtf=genome_gtf,
    output:
        "resources/rrna_{ref}.bed",
    log:
        "resources/logs/rrna_{ref}.bed.rRNA_bed.log",
    cache: True,
    shell:
        """
        python3 scripts/gtf_to_rrna_bed.py -g {input.gtf} -o {output} > {log} 2>&1
        """

