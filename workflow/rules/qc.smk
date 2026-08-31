
rule targetcov_bed:
    input:
        bam= "{runid}/results/reads/star/mrkdup/{sample}.bam",
        bed= config["bed"].split('.')[0] + ".srtd.bed",     #"resources/twist_rna_exome_target_regions_hg38_annotated.srtd.bed",
        genord="resources/genome.txt"
    output:
        "{runid}/results/qc/bedtools/{sample}/bcov.tsv"
    conda:
      "../envs/hts.yaml"
    log: "{runid}/logs/bedtools/{sample}.cov.log"
    threads: 57
    resources:
        mem_gb = 242
    priority: 2
    shell:"""
        bedtools coverage -a {input.bed} -b {input.bam} -g {input.genord} -sorted -d > {output} 2> {log}
        #coverageBed -b {input.bam} -a {input.bed} -split -d -sorted > {output} 2> {log}
        """

rule cov_n10:
    input:
        cov = "{runid}/results/qc/bedtools/{sample}/bcov.tsv",
        genes = "resources/winters_and_cegat_genes.tsv"
    output:
        "{runid}/results/qc/n10_cov/{sample}.n10.tsv"
    conda:
      "../envs/hts.yaml"
    log: "{runid}/logs/bedtools/{sample}.cov_n10.log"
    threads: 6
    shell: """
        python scripts/clinFuseCov.py -c {input.cov} -d {input.genes} > {output} 2> {log}
        """

rule targetcov_perc:
    input:
        "{runid}/results/qc/bedtools/{sample}/bcov.tsv"
    output:
        "{runid}/results/qc/bedtools/{sample}/bcov.perc.tsv"
    log: "{runid}/logs/bedtools/{sample}.awk.log"
    conda:
      "../envs/hts.yaml"
    threads: 14
    priority: 0
    shell: """
        awk '{{if ($6 > 0) sum += 1}} END {{print (sum / NR) * 100}}' {input} | tr ',' '.' > {output} 2> {log}
        """


rule clumpify_opt_dup:
    input:
        r1= "{runid}/results/bcl2fq/{sample}.R1_001.fastq.gz",
        r2= "{runid}/results/bcl2fq/{sample}.R2_001.fastq.gz"
    output:
        r1= temp("{runid}/results/clumpify/opt/{sample}_R1.fq.gz"),
        r2= temp("{runid}/results/clumpify/opt/{sample}_R2.fq.gz"),
        stats= touch("{runid}/results/clumpify/opt/{sample}.clump.stats.txt")
    wildcard_constraints:
        sample= common_constraint
    log:
        "{runid}/logs/clumpify/opt/{sample}.log"
    params:
        extra= "optical",
        dupedist= 12000,   # for NovaSeq. If using NextSeq, then use 40
        subs= 2
    resources:
        mem_gb= 31
    threads: 12
    shell:"""
        clumpify.sh Xmx31g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} dupedist={params.dupedist} dedupe {params.extra} k=23 passes=2 subs={params.subs}  2> {log}
        """

rule clumpify_pcr_dup:
    input:
        r1= "{runid}/results/bcl2fq/{sample}.R1_001.fastq.gz",
        r2= "{runid}/results/bcl2fq/{sample}.R2_001.fastq.gz"
    output:
        r1= temp("{runid}/results/clumpify/pcr/{sample}_R1.fq.gz"),
        r2= temp("{runid}/results/clumpify/pcr/{sample}_R2.fq.gz"),
        stats= touch("{runid}/results/clumpify/pcr/{sample}.clump.stats.txt")
    wildcard_constraints:
        sample= common_constraint
    log:
        "{runid}/logs/clumpify/pcr/{sample}.log"
    params:
        extra= "",
        subs= 2
    resources:
        mem_gb= 31
    threads: 12
    shell:"""
        clumpify.sh Xmx31g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} dedupe {params.extra} k=23 passes=2 subs={params.subs}  2> {log}
        """


rule fastqc_bbmerged:
    input:
        "{runid}/results/bbmerge/{sample}.fq.gz"
    output:
        html="{runid}/results/qc/bbmerge/fastqc/{sample}.html",
        zip="{runid}/results/qc/bbmerge/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    conda:
      "../envs/hts.yaml"
    params:
        extra = "--quiet"
    log:
        "{runid}/logs/qc/fastqc_bbmerge/{sample}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.1.1/bio/fastqc"


rule multiqc_bbmerged:
    input:
        expand("{runid}/results/qc/bbmerge/fastqc/{sample}.html", runid= runid, sample= idkeys)
    output:
        "{runid}/results/qc/fastqc_multiqc_report.html"
    conda:
      "../envs/hts.yaml"
    params:
        extra="--zip-data-dir"
    log:
        "{runid}/logs/qc/multiqc_bbmerged.log"
    wrapper:
        "v2.12.0/bio/multiqc"



rule fastqc_cons:
    input:
        "{runid}/results/reads/consensus/{sample}.{read}.fq"
    output:
        html="{runid}/results/reads/consensus/qc/fastqc/{sample}.{read}.html",
        zip="{runid}/results/reads/consensus/qc/fastqc/{sample}.{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    conda:
      "../envs/hts.yaml"
    params:
        extra = "--quiet"
    log:
        "{runid}/logs/qc/fastqc_cons/{sample}.{read}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.1.1/bio/fastqc"




rule multiqc_cons:
    input:
        expand("{runid}/results/reads/consensus/qc/fastqc/{sample}.{read}.html", runid= runid, sample= idkeys, read= [1,2])
    output:
        "{runid}/results/reads/consensus/qc/multiqc_report.html"
    conda:
      "../envs/hts.yaml"
    params:
        extra="--zip-data-dir"
    log:
        "{runid}/logs/qc/multiqc_umi.log"
    wrapper:
        "v2.1.1/bio/multiqc"







## RSEQC


rule rseqc_gtf2bed:
    input:
        "resources/genome.gtf",
    output:
        bed="{runid}/results/qc/rseqc/annotation.bed",
        db="{runid}/results/qc/rseqc/annotation.db",
    conda:
      "../envs/hts.yaml"
    log:
        "{runid}/logs/rseqc/rseqc_gtf2bed.log",
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="{runid}/results/reads/star/{sample}.bam",
        bai="{runid}/results/reads/star/{sample}.bam.bai",
        bed="resources/annotation.bed", #"{runid}/results/qc/rseqc/annotation.bed",
    output:
        "{runid}/results/qc/rseqc/{sample}.junctionanno.junction.bed",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_junction_annotation/{sample}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="{runid}/results/reads/star/{sample}.bam",
        bed="resources/annotation.bed",  #"{runid}/results/qc/rseqc/annotation.bed",
    output:
        "{runid}/results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1"


rule rseqc_stat:
    input:
        bam="{runid}/results/reads/star/{sample}.bam",
    output:
        "{runid}/results/qc/rseqc/{sample}.stats.txt",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_stat/{sample}.log",
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="{runid}/results/reads/star/{sample}.bam",
        bed="resources/annotation.bed",   #"{runid}/results/qc/rseqc/annotation.bed",
    output:
        "{runid}/results/qc/rseqc/{sample}.infer_experiment.txt",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_infer/{sample}.log",
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="{runid}/results/reads/star/{sample}.bam",
        bed="resources/annotation.bed",   #"{runid}/results/qc/rseqc/annotation.bed",
    output:
        "{runid}/results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_innerdis/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="{runid}/results/reads/star/{sample}.bam",
        bed="resources/annotation.bed",     #"{runid}/results/qc/rseqc/annotation.bed",
    output:
        "{runid}/results/qc/rseqc/{sample}.readdistribution.txt",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_readdis/{sample}.log",
    #conda:
    #     "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "{runid}/results/reads/star/{sample}.bam",
    output:
        "{runid}/results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_readdup/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "{runid}/results/reads/star/{sample}.bam",
    output:
        "{runid}/results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
    conda:
      "../envs/hts.yaml"
    priority: 1
    log:
        "{runid}/logs/rseqc/rseqc_readgc/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqcRSeQC:
    input:
        expand(
            "{runid}/results/reads/star/{sample}.bam",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{sample}.junctionanno.junction.bed",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{sample}.stats.txt",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{sample}.readdistribution.txt",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
            sample= wts_samples, runid= runid,
        ),
        expand(
            "{runid}/logs/rseqc/rseqc_junction_annotation/{sample}.log",
            sample= wts_samples, runid= runid,
        ),
    output:
        "{runid}/results/qc/rseqc_multiqc_report.html",
    conda:
      "../envs/hts.yaml"
    params:
        extra= "",
        #use_input_files_only= True,
    log:
        "{runid}/logs/rseqc_multiqc.log",
    wrapper:
        "v1.23.1/bio/multiqc"


rule deeptools_multiBamSummary:
    input:
        bams = expand("{runid}/results/reads/star/mrkdup/{sample}.bam", runid= runid, sample= idkeys),
        bais = expand("{runid}/results/reads/star/mrkdup/{sample}.bam.bai", runid= runid, sample= idkeys),
        bed = "resources/genome.gtf"
    output:
        "{runid}/results/qc/deeptools/bamSummary.npz"
    conda:
      "../envs/hts.yaml"
    params:
        binSize = 100,
        extra = "--ignoreDuplicates --centerReads --smartLabels"    # see https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html
    log: "{runid}/logs/deeptools/multiBamSummary.log"
    threads: 24
    shell: """
        multiBamSummary BED-file --BED {input.bed} -p {threads} {params.extra} --bamfiles {input.bams} -o {output} 2> {log}
        """
        # here, you can use BED-file instead of bins mode (then, you need to supply the --BED {input.bed}
        # multiBamSummary bins -bs {params.binSize} -p {threads} {params.extra} --bamfiles {input} -o {output} 2> {log}

rule deeptools_plotCor:
    input:
      "{runid}/results/qc/deeptools/bamSummary.npz"
    output:
      "{runid}/results/qc/deeptools/cor_plot.png"
    conda:
      "../envs/hts.yaml"
    params:
      extra = "-p heatmap -c pearson"   # spearman or pearson (i.e. more robust or more sensible?) see https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html
    log: "{runid}/logs/deeptools/plotCor.log"
    shell: """
      plotCorrelation -in {input} {params.extra} -o {output} 2> {log}
      """

rule deeptools_plotPCA:
    input:
      "{runid}/results/qc/deeptools/bamSummary.npz"
    output:
      "{runid}/results/qc/deeptools/pca_plot.png"
    conda:
      "../envs/hts.yaml"
    params:
      extra = "--transpose"    # see https://deeptools.readthedocs.io/en/develop/content/tools/plotPCA.html
    log: "{runid}/logs/deeptools/plotPCA.log"
    shell: """
      plotPCA -in {input} {params.extra} -o {output} 2> {log}
      """
