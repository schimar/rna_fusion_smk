rule rsync:
    input:
        fastqc = expand("{runid}/results/qc/fastqc_multiqc_report.html", runid= runid), 
        rseqc = expand("{runid}/results/qc/rseqc_multiqc_report.html", runid= runid),
        arriba = expand("{runid}/results/arriba/{sample}/fusions.clinout.fltrd.ex.tsv", runid= runid, sample= idkeys),#, unit=units.itertuples()),  #, sample= idkeys[0]),
        rmats = expand("{runid}/results/rmats/{sample}.exonSkip.tsv", runid= runid, sample= idkeys), #unit=units.itertuples()), #, sample= idkeys[0]), 
        bcov = expand("{runid}/results/qc/bedtools/{sample}/bcov.tsv", runid= runid, sample= idkeys),
        n10cov = expand("{runid}/results/qc/n10_cov/{sample}.n10.tsv", runid= runid, sample= idkeys),
    output:
        touch("{bcldir}{analysis_path}/run.done")
    #log:
    #    "{runid}/logs/rsync.log"
    wildcard_constraints:
        sample = common_constraint,
        analysis_path = "analysis_rna.*"
    params:
        bcldir = config["bcldir"],
        runid = config['runid'],
        final_dest = config['bcldir'] + config['analysis_path'],
        log = config['bcldir'] + "{analysis_path}/logs/rsync.log"
    shell:"""
        rsync -rDvz {params.runid}/* {params.final_dest}/ --log-file={params.log}
        """

 ##dt=$(date '+%d%m%Y_%H%M')
        ##mkdir -p {bcldir}/analysis_rna_${{dt}}/qc/{{fastqc,rseqc}} {bcldir}/analysis_rna_${{dt}}/{{arriba,rmats}} &&

  #rsync -rDvz {{input.rseqc} {bcldir}/analysis_rna_${{dt}}/qc/rseqc/ &&
        #rsync -rDvz {input.fastqc} {bcldir}/analysis_rna_${{dt}}/qc/fastqc/ &&
        #rsync -rDvz {runid}/logs/{log} {bcldir}/analysis_rna_${{dt}}/
        #rsync -rDvz {input.arriba} {bcldir}/analysis_rna_${{dt}}/



