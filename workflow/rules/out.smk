rule rsync:
    input:
        fastqc = expand("{runid}/results/qc/fastqc_multiqc_report.html", runid= runid), 
        rseqc = expand("{runid}/results/qc/rseqc_multiqc_report.html", runid= runid),
        arriba = expand("{runid}/results/arriba/{sample}/fusions.clinout.fltrd.ex.tsv", runid= runid, sample= idkeys),#, unit=units.itertuples()),  #, sample= idkeys[0]),
        rmats = expand("{runid}/results/rmats/{sample}.exonSkip.tsv", runid= runid, sample= idkeys), #unit=units.itertuples()), #, sample= idkeys[0]), 
        bcov = expand("{runid}/results/qc/bedtools/{sample}/bcov.tsv", runid= runid, sample= idkeys),
        n10cov = expand("{runid}/results/qc/n10_cov/{sample}.n10.tsv", runid= runid, sample= idkeys),
    output:
        touch("{bcldir}analysis_rna/run.done")
    #log:
    #    "{runid}/logs/rsync.log"
    wildcard_constraints:
        sample = common_constraint
    params:
        bcldir = config["bcldir"],
        runid = config['runid'],
        log = config['bcldir'] + "analysis_rna/logs/rsync.log"
    shell:"""
        rsync -rDvz {params.runid}/* {params.bcldir}analysis_rna/ --log-file={params.log}
        """

 ##dt=$(date '+%d%m%Y_%H%M')
        ##mkdir -p {bcldir}/analysis_rna_${{dt}}/qc/{{fastqc,rseqc}} {bcldir}/analysis_rna_${{dt}}/{{arriba,rmats}} &&

  #rsync -rDvz {{input.rseqc} {bcldir}/analysis_rna_${{dt}}/qc/rseqc/ &&
        #rsync -rDvz {input.fastqc} {bcldir}/analysis_rna_${{dt}}/qc/fastqc/ &&
        #rsync -rDvz {runid}/logs/{log} {bcldir}/analysis_rna_${{dt}}/
        #rsync -rDvz {input.arriba} {bcldir}/analysis_rna_${{dt}}/



