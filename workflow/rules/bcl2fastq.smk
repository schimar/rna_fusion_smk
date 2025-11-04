#ruleorder: bcl2fq > cat_fq1 
#ruleorder: bcl2fq > cat_fq2
rule bcl_convert:
    input:
      bcldir= config['bcldir'],
    output:
      "{runid}/results/bcl2fq/Logs/FastqComplete.txt",
    params:
      outdir= "{runid}/results/bcl2fq/",
      extra= ""
    log:
      "{runid}/logs/bcl_convert.log",
    threads: 56
    shell:"""
      nohup bcl-convert -f --bcl-input-directory {input.bcldir} --output-directory {params.outdir} --sample-sheet {input.bcldir}/SampleSheet_rna.csv --bcl-sampleproject-subdirectories false --bcl-num-parallel-tiles 1 --bcl-num-conversion-threads {threads} > {log} 2>&1
      """
# --shared-thread-odirect-output true 

rule bcl2fq:
    input:
        bcldir= config['bcldir'],
        #bcldir = {bcldir}
        #sampleSheet = "{bcldir}/SampleSheet_rna.csv",
    output:
        #expand("{runid}/results/bcl2fq/{sample}_{read}_001.fastq.gz", runid= runid, sample = idkeys, lane = ['L001', 'L002'], read = ['R1', 'R2']), #, runid = config['runID']),
        "{runid}/results/bcl2fq/Reports/html/tree.html",
#    wildcard_constraints:
#        runid = "[0-9A-Za-z\/]+[^routine]"
    conda:
      "../envs/bcl2fastq.yaml"
    params:
        outdir = "{runid}/results/bcl2fq/",
        units = "{runid}/units.tsv"
    log:
        "{runid}/logs/bcl2fastq.log",
    threads: 32
    shell:
        """
        nohup bcl2fastq --runfolder-dir {input.bcldir} --output-dir {params.outdir} --sample-sheet {input.bcldir}/SampleSheet_rna.csv -p 32 -r 4 -w 4 > {log} 2>&1
        #rm {params.units}
        """
 #{runid}/logs/bcl2fastq.log 2>&1


## expand("../fq/{{runid}}_S0_L001_{readid}_001.fastq.gz", readid=config['readids'])
## I will probably have to use the expand in the output, so I can name the files in rule all! 
   
#   "{sample}_{read}_001.fastq.gz"
rule cat_lanes:
    input:
        fq = expand(config["runid"] + "/results/bcl2fq/{{sample}}_{lane}_R{{read}}_001.fastq.gz", lane=config["lanes"])
        # with 1 lane 
        #fq = expand(config["runid"] + "/results/bcl2fq/{{sample}}_L001_R{{read}}_001.fastq.gz")
    output:
        fq = temp("{runid}/results/reads/cat/{sample}_R{read}.fq.gz")
    conda:
      "../envs/hts.yaml"
    threads: 1
    resources:
        #mem_mb=100,
    #benchmark:
    #    "{runid}/benchmarks/cat_lanes/{sample}_{read}.tsv"
    log:
        "{runid}/logs/cat_lanes/{sample}_{read}.log"
    shell:
        """
        cat {input.fq} > {output.fq} 2> {log}
        """


rule bbmerge_fqs:
    input:
        fq1 = "{runid}/results/reads/cat/{sample}_R1.fq.gz",
        fq2 = "{runid}/results/reads/cat/{sample}_R2.fq.gz",
        #fq1 = "{runid}/results/bcl2fq/{sample}_R1_001.fastq.gz",
        #fq2 = "{runid}/results/bcl2fq/{sample}_R2_001.fastq.gz",        #fq1 = "/mnt/sda/rnaSeq/runs/231025/results/bcl2fq/cat/{sample}_R1.fq.gz",
        #fq2 = "/mnt/sda/rnaSeq/runs/231025/results/bcl2fq/cat/{sample}_R2.fq.gz",
    output:
        out = temp("{runid}/results/bbmerge/{sample}.fq.gz"),
        outu1 = temp("{runid}/results/bbmerge/outu/{sample}.1.fq.gz"),
        outu2 = temp("{runid}/results/bbmerge/outu/{sample}.2.fq.gz"),
        hist = "{runid}/results/bbmerge/{sample}.hist.txt"
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log: "{runid}/logs/bbmerge/{sample}.log"
    priority: 1
    resources:
        mem_gb=60
    threads: 12
    shell:"""
        bbmerge-auto.sh -Xmx60g in1={input.fq1} in2={input.fq2} outm={output.out} outu1={output.outu1} outu2={output.outu2} ihist={output.hist} ecct extend2=20 iterations=5 > {log} 2>&1
        """








#rule cat_fq1:
#    input:
#        "{runid}/results/bcl2fq/Reports/html/tree.html"
#    output: 
#        fq = temp("{runid}/results/bcl2fq/cat/{sample}_R1.fq.gz"),
#        #fq = "/mnt/illumina/230627_A01358_0051_BHGCJ5DRX3/Fastq/{runid}/cat/{sample}_R1.fq.gz",
#        done = touch("{runid}/results/bcl2fq/cat/{sample}.cat.1.done")
#    log: "{runid}/logs/cat_fq/{sample}_1.log"
#    shell: """
#        cat {runid}/results/bcl2fq/{wildcards.sample}_L*_R1_001.fastq.gz > {output.fq} 2> {log}
#        """
#
#rule cat_fq2:
#    input:
#        "{runid}/results/bcl2fq/Reports/html/tree.html"
#    output: 
#        fq = temp("{runid}/results/bcl2fq/cat/{sample}_R2.fq.gz"),
#        #fq = "/mnt/illumina/230627_A01358_0051_BHGCJ5DRX3/Fastq/{runid}/cat/{sample}_R2.fq.gz",
#        done = touch("{runid}/results/bcl2fq/cat/{sample}.cat.2.done")
#    log: "{runid}/logs/cat_fq/{sample}_2.log"
#    shell: """
#        cat {runid}/results/bcl2fq/{wildcards.sample}_L*_R2_001.fastq.gz > {output.fq} 2> {log}
#        """


