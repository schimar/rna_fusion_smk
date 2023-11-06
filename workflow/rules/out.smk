rule mv_results:
    input:
        rseqc = expand("{runid}/results/bbmerge/qc/multiqc_report.html", runid= runid), 
        fastqc = expand("{runid}/results/qc/multiqc_report.html", runid= runid),
        arriba = expand("{runid}/results/arriba/{sample}/fusions.clinout.fltrd.ex.tsv", runid= runid, sample= idkeys),#, unit=units.itertuples()),  #, sample= idkeys[0]),
        rmats = expand("{runid}/results/rmats/{sample}/{sample}.exonSkip.tsv", runid= runid, sample= idkeys), #unit=units.itertuples()), #, sample= idkeys[0]), 
    output:
        touch("{runid}/run.done")
    wildcard_constraints:
        sample = common_constraint
    #params:
        #bcldir = "{bcldir}"
    shell:"""
        dt=$(date '+%d%m%Y_%H%M')
        mkdir -p {bcldir}/analysis_rna_${{dt}}/qc/{{fastqc,rseqc}} {bcldir}/analysis_rna_${{dt}}/{{arriba,rmats}} &&
        rsync -rDvz {input.rseqc} {bcldir}/analysis_rna_${{dt}}/qc/rseqc/ &&
        rsync -rDvz {input.fastqc} {bcldir}/analysis_rna_${{dt}}/qc/fastqc/ &&
        rsync -rDvz {runid}/logs/{log} {bcldir}/analysis_rna_${{dt}}/
        #rsync -rDvz {input.arriba} {bcldir}/analysis_rna_${{dt}}/
        """




#rule arriba:
#    input:
#        bam="{runid}/results/reads/star/mrkdup/{sample}.bam",
#        genome="resources/genome.fa",
#        annotation="resources/genome.gtf",
#        # optional: # A custom tsv containing identified artifacts, such as read-through fusions of neighbouring genes.
#        #blacklist="resources/blacklist_hg38_GRCh38_v2.3.0.tsv.gz",
#        #known_fusions="resources/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz",
#        # default blacklists are selected via blacklist parameter
#        # see https://arriba.readthedocs.io/en/latest/input-files/#blacklist
#        #custom_blacklist=[],
#    output:
#        fusions="{runid}/results/arriba/{sample}/fusions.tsv",
#        discarded="{runid}/results/arriba/{sample}/fusions.discarded.tsv",
#        #done="{runid}/results/arriba/{sample}/arriba.done",
#    params:
#        # required if blacklist or known_fusions is set
#        genome_build="GRCh38",
#        default_blacklist=True,
#        default_known_fusions=True,
#        extra="-u",   #alignIntronMax",
#    log:
#        "{runid}/logs/arriba/{sample}.log",
#    threads: 1
#    wrapper:
#        "v1.23.4/bio/arriba"


#rule get_clinFuse:
#    input:
#        fus = "{runid}/results/arriba/{sample}/fusions.tsv",
#        clinTab = "resources/s2_winters2018.tsv",
#    output:
#        clinout = "{runid}/results/arriba/{sample}/fusions.clinout.tsv",
#        nonclinout = "{runid}/results/arriba/{sample}/fusions.nonclinout.tsv",
#    shell:"""
#        scripts/clinFuse.py -f {input.fus} -c {input.clinTab}
#        """
#
#
#rule filter_arriba:
#    input:
#        "{runid}/results/arriba/{sample}/fusions.{clin_nonclin}.tsv",
#    output:
#        "{runid}/results/arriba/{sample}/fusions.{clin_nonclin}.fltrd.tsv",
#    #wildcard_constraints:
#    #    sample = common_constraint
#    log:
#        "{runid}/logs/fltr_arriba/{sample}/{clin_nonclin}.log"
#    shell:"""
#        scripts/fusionfltr.py -i {input} > {output} 2> {log}
#        """
#
#
#rule getExons_arriba:
#    input:
#        tsv = "{runid}/results/arriba/{sample}/fusions.clinout.fltrd.tsv",
#        h5 = "resources/exons23.h5",
#    output:
#        "{runid}/results/arriba/{sample}/fusions.clinout.fltrd.ex.tsv",
#    #wildcard_constraints:
#    #    sample = common_constraint
#    log:
#        "{runid}/logs/getExons_arriba/{sample}.log"
#    shell:"""
#        scripts/fusGetExons.py -i {input.tsv} -d {input.h5} > {output} 2> {log}
#        """


