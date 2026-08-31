
rule arriba:
    input:
        bam="{runid}/results/bam/{sample}.bam",
        genome=genome_fa,
        annotation=genome_gtf,
        # optional: # A custom tsv containing identified artifacts, such as read-through fusions of neighbouring genes.
        #blacklist="resources/blacklist_hg38_GRCh38_v2.3.0.tsv.gz",
        #known_fusions="resources/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz",
        # default blacklists are selected via blacklist parameter
        # see https://arriba.readthedocs.io/en/latest/input-files/#blacklist
        custom_blacklist = arriba_blacklist,  # set in config.yaml; #[],
    output:
        fusions="{runid}/results/fusions/{sample}.fusions.tsv",
        discarded="{runid}/results/fusions/{sample}.fusions.discarded.tsv",
        #done="{runid}/results/fusions/{sample}.arriba.done",
    params:
        # required if blacklist or known_fusions is set
        genome_build="GRCh38",
        default_blacklist=False,
        default_known_fusions=True,
        extra="-u",   #alignIntronMax",
    log:
        "{runid}/results/fusions/{sample}.arriba.log",
    threads: 1
    wrapper:
        "v1.23.4/bio/arriba"


rule get_clinFuse:
    input:
        fus = "{runid}/results/arriba/{sample}/fusions.tsv",
        clinTab = "resources/winters_and_cegat_genes.tsv",
    output:
        clinout = "{runid}/results/arriba/{sample}/fusions.clinout.tsv",
        nonclinout = "{runid}/results/arriba/{sample}/fusions.nonclinout.tsv",
    conda:
      "../envs/hts.yaml"
    log: 
        "{runid}/logs/clinFuse/{sample}.log"
    shell:"""
        python scripts/clinFuse.py -f {input.fus} -c {input.clinTab} 2> {log}
        """


rule filter_arriba:
    input:
        fus = "{runid}/results/arriba/{sample}/fusions.{clin_nonclin}.tsv",
        bl = "resources/bl_exonSkip.tsv",
    output:
        "{runid}/results/arriba/{sample}/fusions.{clin_nonclin}.fltrd.tsv",
    conda:
      "../envs/hts.yaml"
    #wildcard_constraints:
    #    sample = common_constraint
    log:
        "{runid}/logs/fltr_arriba/{sample}/{clin_nonclin}.log"
    shell:"""
        python scripts/fusionfltr.py -i {input.fus} -b {input.bl} > {output} 2> {log}
        """

rule gtf2exon_h5:
    input:
      gtf = "resources/genome.gtf"
    output:
      h5 = "resources/exons_from_gtf.h5"
    conda:
      "../envs/hts.yaml"
    log:
      "logs/gtf2exon_h5.log"
    shell: """
      python scripts/gtf2exon_h5.py -g {input.gtf} -h5 {output.h5}
      """


rule getExons_arriba:
    input:
        tsv = "{runid}/results/arriba/{sample}/fusions.clinout.fltrd.tsv",
        h5 = "resources/exons_from_gtf.h5",
    output:
        "{runid}/results/arriba/{sample}/fusions.clinout.fltrd.ex.tsv",
    conda:
      "../envs/hts.yaml"
    #wildcard_constraints:
    #    sample = common_constraint
    log:
        "{runid}/logs/getExons_arriba/{sample}.log"
    shell:"""
        python scripts/fus_get_exonNo.py -i {input.tsv} -d {input.h5} > {output} 2> {log}
        """



