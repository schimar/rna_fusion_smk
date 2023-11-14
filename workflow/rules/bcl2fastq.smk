#ruleorder: bcl2fq > cat_fq1 
#ruleorder: bcl2fq > cat_fq2


rule bcl2fq:
    input:
        bcldir = {bcldir} #"/mnt/illumina/231025_A01272_0063_AHK5CGDRX3/",
        #bcldir = "/mnt/routine/230822_A01358_0062_AHJWH7DRX3/",
        #bcldir = "/mnt/routine/230915_A01358_0066_AHJWHKDRX3/",
        #bcldir = "/mnt/illumina/231011_A01272_0061_AHK7N7DRX3/",
        #230426_A01272_0045_AH5CT5DRX3/",
    output:
        #expand("{runid}/results/bcl2fq/{sample}_{read}_001.fastq.gz", runid= runid, sample = idkeys, lane = ['L001', 'L002'], read = ['R1', 'R2']), #, runid = config['runID']),
        "{runid}/results/bcl2fq/Reports/html/tree.html",
#    wildcard_constraints:
#        runid = "[0-9A-Za-z\/]+[^routine]"
    params:
        outdir = "{runid}/results/bcl2fq/",
    #log:
    #    "{runid}/logs/bcl2fastq.log",
    threads: 32
    shell:
        """
        nohup bcl2fastq --runfolder-dir {input.bcldir} --output-dir {params.outdir} -p 18 -r 3 -w 3 > {runid}/logs/bcl2fastq.log 2>&1
        """



## expand("../fq/{{runid}}_S0_L001_{readid}_001.fastq.gz", readid=config['readids'])
## I will probably have to use the expand in the output, so I can name the files in rule all! 
   
#   "{sample}_{read}_001.fastq.gz"

rule cat_fq1:
    input:
        "{runid}/results/bcl2fq/Reports/html/tree.html"
    output: 
        fq = temp("{runid}/results/bcl2fq/cat/{sample}_R1.fq.gz"),
        #fq = "/mnt/illumina/230627_A01358_0051_BHGCJ5DRX3/Fastq/{runid}/cat/{sample}_R1.fq.gz",
        done = touch("{runid}/results/bcl2fq/cat/{sample}.cat.1.done")
    log: "{runid}/logs/cat_fq/{sample}_1.log"
    shell: """
        cat {runid}/results/bcl2fq/{wildcards.sample}_L*_R1_001.fastq.gz > {output.fq} 2> {log}
        """

rule cat_fq2:
    input:
        "{runid}/results/bcl2fq/Reports/html/tree.html"
    output: 
        fq = temp("{runid}/results/bcl2fq/cat/{sample}_R2.fq.gz"),
        #fq = "/mnt/illumina/230627_A01358_0051_BHGCJ5DRX3/Fastq/{runid}/cat/{sample}_R2.fq.gz",
        done = touch("{runid}/results/bcl2fq/cat/{sample}.cat.2.done")
    log: "{runid}/logs/cat_fq/{sample}_2.log"
    shell: """
        cat {runid}/results/bcl2fq/{wildcards.sample}_L*_R2_001.fastq.gz > {output.fq} 2> {log}
        """


rule bbmerge_fqs:
    input:
        fq1 = "{runid}/results/bcl2fq/cat/{sample}_R1.fq.gz",
        fq2 = "{runid}/results/bcl2fq/cat/{sample}_R2.fq.gz",
        #fq1 = "/mnt/sda/rnaSeq/runs/231025/results/bcl2fq/cat/{sample}_R1.fq.gz",
        #fq2 = "/mnt/sda/rnaSeq/runs/231025/results/bcl2fq/cat/{sample}_R2.fq.gz",
    output:
        out = temp("{runid}/results/bbmerge/{sample}.fq.gz"),
        outu1 = temp("{runid}/results/bbmerge/outu/{sample}.1.fq.gz"),
        outu2 = temp("{runid}/results/bbmerge/outu/{sample}.2.fq.gz"),
        hist = "{runid}/results/bbmerge/{sample}.hist.txt"
    wildcard_constraints:
        sample = common_constraint
    log: "{runid}/logs/bbmerge/{sample}.log"
    resources:
        mem_gb=46
    threads: 8
    shell:"""
        bbmerge-auto.sh -Xmx46g in1={input.fq1} in2={input.fq2} outm={output.out} outu1={output.outu1} outu2={output.outu2} ihist={output.hist} ecct extend2=20 iterations=5 > {log} 2>&1
        """

