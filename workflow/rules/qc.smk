## bb 


# bbmerge.sh in1=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R1.fq.gz in2=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R2.fq.gz out=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/H22032902-25ng-rep1_S13.fq.gz outu1=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/outu/H22032902-25ng-rep1_S13_1.fq.gz outu2=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/outu/H22032902-25ng-rep1_S13_2.fq.gz       
    
# bbmerge-auto.sh in1=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R1.fq.gz in2=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R2.fq.gz out=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/kmer/H22032902-25ng-rep1_S13.fq.gz outu1=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/kmer/outu/H22032902-25ng-rep1_S13_1.fq.gz outu2=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/kmer/outu/H22032902-25ng-rep1_S13_2.fq.gz ihist=ihist_auto.txt ecct extend2=20 iterations=5

rule fastqc_bbmerged:
    input:
        "{runid}/results/bbmerge/{sample}.fq.gz"
    output:
        html="{runid}/results/qc/bbmerge/fastqc/{sample}.html",
        zip="{runid}/results/qc/bbmerge/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
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
            "{runid}/results/reads/star/{unit.sample_name}.bam",
            unit=units.itertuples(), runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.junctionanno.junction.bed",
            unit=units.itertuples(), runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.junctionsat.junctionSaturation_plot.pdf",
            unit=units.itertuples(), runid= runid,
        ),
        #expand(
        #    "{runid}/results/qc/rseqc/{unit.sample_name}.infer_experiment.txt",
        #    unit=units.itertuples(), runid= runid,
        #),
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.stats.txt",
            unit=units.itertuples(), runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.inner_distance_freq.inner_distance.txt",
            unit=units.itertuples(), runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.readdistribution.txt",
            unit=units.itertuples(), runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.readdup.DupRate_plot.pdf",
            unit=units.itertuples(), runid= runid,
        ),
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.readgc.GC_plot.pdf",
            unit=units.itertuples(), runid= runid,
        ),
        expand(
            "{runid}/logs/rseqc/rseqc_junction_annotation/{unit.sample_name}.log",
            unit=units.itertuples(), runid= runid,
        ),
    output:
        "{runid}/results/qc/rseqc_multiqc_report.html",
    params:
        extra= "",
        #use_input_files_only= True,
    log:
        "{runid}/logs/rseqc_multiqc.log",
    wrapper:
        "v1.23.1/bio/multiqc"


