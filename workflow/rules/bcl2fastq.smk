rule bcl2fq:
    input:
        #rundir = "/mnt/illumina/230209_A01272_0035_BHTVFGDRX2/",
        rundir = "/mnt/illumina/230627_A01358_0051_BHGCJ5DRX3/",
        #230426_A01272_0045_AH5CT5DRX3/",
    output:
        #expand("{runid}/results/bcl2fq/{sample}_{read}_001.fastq.gz", runid= runid, sample = idkeys, lane = ['L001', 'L002'], read = ['R1', 'R2']), #, runid = config['runID']),
        "{runid}/results/bcl2fq/Reports/html/tree.html",
    params:
        outdir = "{runid}/results/bcl2fq/",
    #log:
    #    "{runid}/logs/bcl2fastq.log",
    threads: 32
    shell:
        """
        nohup bcl2fastq --runfolder-dir {input.rundir} --output-dir {params.outdir} -p 18 -r 3 -w 3 > {runid}/logs/bcl2fastq.log 2>&1
        """



## expand("../fq/{{runid}}_S0_L001_{readid}_001.fastq.gz", readid=config['readids'])
## I will probably have to use the expand in the output, so I can name the files in rule all! 
    
#   "{sample}_{read}_001.fastq.gz"

rule cat_fq1:
    output: 
        fq = "{runid}/results/bcl2fq/cat/{sample}_R1.fq.gz",
        #fq = "/mnt/illumina/230627_A01358_0051_BHGCJ5DRX3/Fastq/{runid}/cat/{sample}_R1.fq.gz",
        done = touch("{runid}/results/bcl2fq/cat/{sample}.cat.1.done")
    log: "{runid}/logs/cat_fq/{sample}_1.log"
    shell: """
        cat {runid}/results/bcl2fq/{wildcards.sample}_L*_R1_001.fastq.gz > {output.fq} 2> {log}
        """

rule cat_fq2:
    output: 
        fq = "{runid}/results/bcl2fq/cat/{sample}_R2.fq.gz",
        #fq = "/mnt/illumina/230627_A01358_0051_BHGCJ5DRX3/Fastq/{runid}/cat/{sample}_R2.fq.gz",
        done = touch("{runid}/results/bcl2fq/cat/{sample}.cat.2.done")
    log: "{runid}/logs/cat_fq/{sample}_2.log"
    shell: """
        cat {runid}/results/bcl2fq/{wildcards.sample}_L*_R2_001.fastq.gz > {output.fq} 2> {log}
        """
