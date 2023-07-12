## bb 

rule bbmerge_fqs:
    input:
        fq1 = "{runid}/results/bcl2fq/cat/{sample}_R1.fq.gz",
        fq2 = "{runid}/results/bcl2fq/cat/{sample}_R2.fq.gz",
    output:
        out = "{runid}/results/bbmerge/{sample}.fq.gz",
        outu1 = "{runid}/results/bbmerge/outu/{sample}.1.fq.gz",
        outu2 = "{runid}/results/bbmerge/outu/{sample}.2.fq.gz",
        hist = "{runid}/results/bbmerge/{sample}.hist.txt"
    wildcard_constraints:
        sample = common_constraint
    log: "{runid}/logs/bbmerge/{sample}.log"
    resources:
        mem_gb=20
    threads: 8
    shell:"""
        bbmerge-auto.sh in1={input.fq1} in2={input.fq2} outm={output.out} outu1={output.outu1} outu2={output.outu2} ihist={output.hist} ecct extend2=20 iterations=5 > {log} 2>&1
        """

# bbmerge.sh in1=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R1.fq.gz in2=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R2.fq.gz out=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/H22032902-25ng-rep1_S13.fq.gz outu1=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/outu/H22032902-25ng-rep1_S13_1.fq.gz outu2=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/outu/H22032902-25ng-rep1_S13_2.fq.gz       
    
# bbmerge-auto.sh in1=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R1.fq.gz in2=/mnt/sda/rnaSeq/runs/230627/results/bcl2fq/cat/H22032902-25ng-rep1_S13_R2.fq.gz out=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/kmer/H22032902-25ng-rep1_S13.fq.gz outu1=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/kmer/outu/H22032902-25ng-rep1_S13_1.fq.gz outu2=/mnt/sda/rnaSeq/runs/230627/results/bbmerge/kmer/outu/H22032902-25ng-rep1_S13_2.fq.gz ihist=ihist_auto.txt ecct extend2=20 iterations=5

rule fastqc_bbmerged:
    input:
        "{runid}/results/bbmerge/{sample}.fq.gz"
    output:
        html="{runid}/results/bbmerge/qc/fastqc/{sample}.html",
        zip="{runid}/results/bbmerge/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "{runid}/logs/qc/fastqc_bbmerge/{sample}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.1.1/bio/fastqc"



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


rule multiqc_bbmerged:
    input:
        expand("{runid}/results/bbmerge/qc/fastqc/{sample}.html", runid= runid, sample= idkeys)
    output:
        "{runid}/results/bbmerge/qc/multiqc_report.html"
    log:
        "{runid}/logs/qc/multiqc_bbmerged.log"
    wrapper:
        "v2.1.1/bio/multiqc"



rule multiqc_cons:
    input:
        expand("{runid}/results/reads/consensus/qc/fastqc/{sample}.{read}.html", runid= runid, sample= idkeys, read= [1,2])
    output:
        "{runid}/results/reads/consensus/qc/multiqc_report.html"
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
        db=temp("{runid}/results/qc/rseqc/annotation.db"),
    log:
        "{runid}/logs/rseqc_gtf2bed.log",
    #conda:
    #    "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="{runid}/results/fq2cons/reads/mapped/{sample}.bam",
        bed="{runid}/results/qc/rseqc/annotation.bed",
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
        bam="{runid}/results/star/{sample}/aligned.out.bam",
        bed="{runid}/results/qc/rseqc/annotation.bed",
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
        "{runid}/results/star/{sample}/aligned.out.bam",
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
        bam="{runid}/results/star/{sample}/aligned.out.bam",
        bed="{runid}/results/qc/rseqc/annotation.bed",
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
        bam="{runid}/results/star/{sample}/aligned.out.bam",
        bed="{runid}/results/qc/rseqc/annotation.bed",
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
        bam="{runid}/results/star/{sample}/aligned.out.bam",
        bed="{runid}/results/qc/rseqc/annotation.bed",
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
        "{runid}/results/star/{sample}/aligned.out.bam",
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
        "{runid}/results/star/{sample}/aligned.out.bam",
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
            "{runid}/results/star/{unit.sample_name}/aligned.out.bam",
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
        expand(
            "{runid}/results/qc/rseqc/{unit.sample_name}.infer_experiment.txt",
            unit=units.itertuples(), runid= runid,
        ),
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
        "{runid}/results/qc/multiqc_report.html",
    log:
        "{runid}/logs/multiqc.log",
    wrapper:
        "v1.23.1/bio/multiqc"


