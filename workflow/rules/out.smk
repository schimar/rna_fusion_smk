rule rsync:
    input:
        multiqc = expand("{runid}/results/quality_control/fastq/multiqc_report.html", runid= runid),
        fusions = expand("{runid}/results/fusions/{sample}.fusions.tsv", runid= runid, sample= wts_samples),
        fusions_disc = expand("{runid}/results/fusions/{sample}.fusions.discarded.tsv", runid= runid, sample= wts_samples),
        splicing = expand("{runid}/results/splicing/{sample}.SE.MATS.JC.txt", runid= runid, sample= wts_samples),
        egfr = expand("{runid}/results/splicing/{sample}.egfr_v3.out", runid= runid, sample= wts_samples),
        clump  = expand("{runid}/results/quality_control/fastq/{sample}.clump.{mode}.stats.txt", runid= runid, sample= wts_samples, mode= ['opt', 'pcr']),
        read_distribution = expand("{runid}/results/quality_control/bam/{sample}.read_distribution.txt", runid=runid, sample=wts_samples),
        globin_read_distribution = expand("{runid}/results/quality_control/bam/{sample}.globin.read_distribution.txt", runid=runid, sample=wts_samples),
        mane_read_distribution = expand("{runid}/results/quality_control/bam/{sample}.mane.readdistribution.txt", runid=runid, sample=wts_samples),
        rrna_contamination = expand("{runid}/results/quality_control/bam/{sample}.rrna_contamination.tsv", runid=runid, sample=wts_samples),
        rrna_custom_content = expand("{runid}/results/quality_control/bam/multiqc_rrna_contamination.json", runid=runid),
        bam_stat = expand("{runid}/results/quality_control/bam/{sample}.bam_stat.txt", runid=runid, sample=wts_samples),
        infer_experiment = expand("{runid}/results/quality_control/bam/{sample}.infer_experiment.txt", runid=runid, sample=wts_samples),
        inner_distance = expand("{runid}/results/quality_control/bam/{sample}.inner_distance.txt", runid=runid, sample=wts_samples),
        junction_annotation = expand("{runid}/results/quality_control/bam/{sample}.junction_annotation.bed", runid=runid, sample=wts_samples),
        junction_saturation = expand("{runid}/results/quality_control/bam/{sample}.junction_saturation.pdf", runid=runid, sample=wts_samples),
        read_duplication = expand("{runid}/results/quality_control/bam/{sample}.read_duplication.pdf", runid=runid, sample=wts_samples),
        read_gc = expand("{runid}/results/quality_control/bam/{sample}.read_gc.pdf", runid=runid, sample=wts_samples),
    output:
        touch("{bcldir}{analysis_path}/run.done"),
    wildcard_constraints:
        analysis_path = "analysis_rna.*",
    params:
        bcldir = config["bcldir"],
        runid = config['runid'],
        final_dest = config['bcldir'] + config['analysis_path'],
        log = config['bcldir'] + config['analysis_path'] + "/rsync.log",
    shell:"""
        mkdir -p {params.final_dest} &&
        rsync -rDvz {params.runid}/results/ {params.final_dest}/ --log-file={params.log}
        """

 ##dt=$(date '+%d%m%Y_%H%M')
        ##mkdir -p {bcldir}/analysis_rna_${{dt}}/qc/{{fastqc,rseqc}} {bcldir}/analysis_rna_${{dt}}/{{arriba,rmats}} &&

  #rsync -rDvz {{input.rseqc} {bcldir}/analysis_rna_${{dt}}/qc/rseqc/ &&
        #rsync -rDvz {input.fastqc} {bcldir}/analysis_rna_${{dt}}/qc/fastqc/ &&
        #rsync -rDvz {runid}/logs/{log} {bcldir}/analysis_rna_${{dt}}/
        #rsync -rDvz {input.arriba} {bcldir}/analysis_rna_${{dt}}/



