rule samtoo_index:
    input:
      "{runid}/results/reads/star/{sample}.bam", 
    output:
      "{runid}/results/reads/star/{sample}.bam.bai",
    shell:
      """
      samtools index {input} 
      """


rule sub_chr7:
    input:
      bam="{runid}/results/reads/star/{sample}.bam", 
      bai="{runid}/results/reads/star/{sample}.bam.bai",
    output:
      chr7="{runid}/results/reads/star/{sample}.chr7.bam",
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
    shell:
      """
      # get the median read length from bam 
      samtools view -F 4 {input.chr7} | awk '{print length($10)}' | sort -u | awk '{ a[i++]=$1; } END { print a[int(i/2)]; }' > {output.medRL}   &&
      # 
      ls {input.chr7} > {output.bamls}
      """


rule rMats:
    input:
        bamls="{runid}/results/reads/star/{sample}/bam.list",
    output:
        directory("{runid}/results/rmats/{sample}/"),
    params:

    log:
        "{runid}/logs/rmats/{sample}.log",
    threads: 12
    shell:  
        """
        rmats.py --b1 {input.bam} --readLength 101 --nthread {threads} --od {output} --tmp 
        --variable-read-length 
        """


rule grepMETse:
    input:
      se="{runid}/results/rmats/{sample}/SE.MATS.JC.txt"
      targets=config['splicingTargets']
    output:
      "{runid}/results/rmats/{sample}/seMET.txt",    
    shell:
        """
        head -1 {input} > {output}" &&
        readarray -t targets < <(cat {input.targets})
        for i in ${targets[@]); do 
          egrep -f {targets} {input} >> {output};
        done
        """

