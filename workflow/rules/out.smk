rule rsync:
    input:
        fastqc = expand("{runid}/results/qc/fastqc_multiqc_report.html", runid= runid),
        arriba = expand("{runid}/results/arriba/{sample}/fusions.tsv", runid= runid, sample= wts_samples),
        arriba_disc = expand("{runid}/results/arriba/{sample}/fusions.discarded.tsv", runid= runid, sample= wts_samples),
        rmats = expand("{runid}/results/rmats/{sample}/SE.MATS.JC.txt", runid= runid, sample= wts_samples),
        egfr = expand("{runid}/results/rmats/{sample}.egfr_v3.out", runid= runid, sample= wts_samples),
        bcov = expand("{runid}/results/qc/bedtools/{sample}/bcov.tsv", runid= runid, sample= wts_samples),
        n10cov = expand("{runid}/results/qc/n10_cov/{sample}.n10.tsv", runid= runid, sample= wts_samples),
        clump  = expand("{runid}/results/clumpify/{mode}/{sample}.clump.stats.txt", runid= runid, sample= wts_samples, mode= ['opt', 'pcr']),
    output:
        touch("{bcldir}{analysis_path}/run.done"),
    wildcard_constraints:
        analysis_path = "analysis_rna.*",
    params:
        bcldir = config["bcldir"],
        runid = config['runid'],
        final_dest = config['bcldir'] + config['analysis_path'],
        log = config['bcldir'] + "{analysis_path}/logs/rsync.log",
    shell:"""
        rsync -rDvz {params.runid}/* {params.final_dest}/ --log-file={params.log}
        """

 ##dt=$(date '+%d%m%Y_%H%M')
        ##mkdir -p {bcldir}/analysis_rna_${{dt}}/qc/{{fastqc,rseqc}} {bcldir}/analysis_rna_${{dt}}/{{arriba,rmats}} &&

  #rsync -rDvz {{input.rseqc} {bcldir}/analysis_rna_${{dt}}/qc/rseqc/ &&
        #rsync -rDvz {input.fastqc} {bcldir}/analysis_rna_${{dt}}/qc/fastqc/ &&
        #rsync -rDvz {runid}/logs/{log} {bcldir}/analysis_rna_${{dt}}/
        #rsync -rDvz {input.arriba} {bcldir}/analysis_rna_${{dt}}/



