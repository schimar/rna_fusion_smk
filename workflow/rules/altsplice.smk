rule samtools_index:
    input:
        "{runid}/results/reads/star/{sample}.bam", 
    output:
      "{runid}/results/reads/star/{sample}.bam.bai",
    log:
      "{runid}/logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"


rule sub_chr7:
    input:
      bam="{runid}/results/reads/star/{sample}.bam", 
      bai="{runid}/results/reads/star/{sample}.bam.bai",
    output:
      chr7="{runid}/results/reads/star/{sample}/chr7.bam",
    wildcard_constraints:
        sample = common_constraint
    log: "{runid}/logs/sub_chr7/{sample}.log"
    shell:
      """
      samtools view -F4 {input.bam} -b chr7 > {output.chr7}
      """


rule rmats_createin:
    input:
        #bam = "{runid}/results/reads/star/{sample}.bam",
        chr7 = "{runid}/results/reads/star/{sample}/chr7.bam",
    output:
        bamls="{runid}/results/reads/star/{sample}/bam.list",
        medRL="{runid}/results/reads/star/{sample}/medianRL.txt",
    wildcard_constraints:
        sample = common_constraint
    log: "{runid}/logs/rmats_createin/{sample}.log"
    shell:
        """
        # get the median read length from bam 
        samtools view -F 4 {input.chr7} | awk '{{print length($10)}}' | sort -u | awk '{{ a[i++]=$1; }} END {{ print a[int(i/2)]; }}' > {output.medRL}   &&
        # 
        ls {input.chr7} > {output.bamls}
        """


rule rMats:
    input:
        bamls = "{runid}/results/reads/star/{sample}/bam.list",
        medRL="{runid}/results/reads/star/{sample}/medianRL.txt",
        gtf = "resources/genome.gtf",
    output:
        direc = directory("{runid}/results/rmats/{sample}/"),
        se = "{runid}/results/rmats/{sample}/SE.MATS.JC.txt",
    wildcard_constraints:
        sample = common_constraint
    params:
        extra = "--variable-read-length --statoff"
    log:
        "{runid}/logs/rmats/{sample}.log",
    threads: 12
    shell:  
        """
        readLen=$(cat {input.medRL} ) &&
        rmats.py --b1 {input.bamls} --readLength ${{readLen}} --nthread {threads} --od {output.direc} --gtf {input.gtf} --tmp {output.direc} {params.extra} > {log} 2>&1
        """


rule grepMETse:
    input:
      se = "{runid}/results/rmats/{sample}/SE.MATS.JC.txt",
      targets = config['splicingTargets']
    output:
      "{runid}/results/rmats/{sample}/seMET.txt",    
    shell:
        """
        head -1 {input.se} > {output} &&
        readarray -t < <(cat {input.targets})
        for i in ${{MAPFILE[@]}}; do 
          egrep -f {input.targets} {input.se} >> {output};
        done
        """

