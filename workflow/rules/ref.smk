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
        "resources/genome.fa",
    output:
        "resources/genome.fa.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.77.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/genome.fa",
    output:
        multiext("resources/genome.fa", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        "logs/bwa-mem2_index/genome.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v1.23.5/bio/bwa-mem2/index"


rule star_index:
    input:
        fasta="resources/genome.fa",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 8
    params:
      extra=lambda wc, input:"--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    wrapper:
        "v1.28.0/bio/star/index"  
        #"0.77.0/bio/star/index"



rule create_dict:
    input:
        "resources/genome.fa",
    output:
        "resources/genome.dict",
    log:
        "logs/picard/create_dict.log",
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
        dict="resources/genome.dict",
    output:
        "resources/{bed}_file_UMI_demo_data_hg38.interval_list",
    log:
        "logs/picard/bedtointervallist/{bed}_file_UMI_demo_data_hg38.log",
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
    config["bed"]
  output:
    "{bed_file_stem}.srtd.bed"
  log: "logs/sort_bed/{bed_file_stem}.log"
  shell:"""
    sort -k1,1V -k2,2n -k3,3n {input} > {output}
    """
