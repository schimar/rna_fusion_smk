ruleorder: samtools_index > egfr_v3
#ruleorder: sub_chr > sambamba_rmdup_sub

rule samtools_index:
    input:
        "{runid}/results/reads/star/mrkdup/{sample}.bam", 
    output:
      "{runid}/results/reads/star/mrkdup/{sample}.bam.bai",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log:
      "{runid}/logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 9  # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"



rule sub_chr:
    input:
      bam = "{runid}/results/reads/star/mrkdup/{sample}.bam", 
      bai = "{runid}/results/reads/star/mrkdup/{sample}.bam.bai",
    output:
      sub = temp("{runid}/results/reads/star/mrkdup/{sample}/{chrom}.dup.bam"),
    conda:
      "../envs/hts.yaml"
    #wildcard_constraints:
    #    sample = common_constraint
    log: "{runid}/logs/sub_{chrom}/{sample}.log",
    shell:
      """
      samtools view -F4 {input.bam} -b {wildcards.chrom} > {output.sub}
      """

rule sambamba_rmdup_sub:
    input:
        "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.dup.bam"
    output:
        "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.rmdup.bam",
    conda:
      "../envs/hts.yaml"
    priority: 20
    params:
        extra="-r"  # optional parameters
    log: "{runid}/logs/sambamba-markdup_subchr/{sample}_{chrom}.log"
    threads: 8
    wrapper:
        "v1.31.1/bio/sambamba/markdup"


rule rmats_createin:
    input:
        bam = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.rmdup.bam",
        #bam = "{runid}/results/reads/star/{sample}/chr7.bam",
    output:
        bamls = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.bam.list",
        medRL = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.medianRL.txt",
    conda:
      "../envs/hts.yaml"
    #wildcard_constraints:
    #    sample = common_constraint
    log: "{runid}/logs/rmats_createin/{sample}.{chrom}.log"
    shell:
        """
        # get the median read length from bam 
        samtools view -F 4 {input.bam} | awk '{{print length($10)}}' | sort -u | awk '{{ a[i++]=$1; }} END {{ print a[int(i/2)]; }}' > {output.medRL}   &&
        # 
        ls {input.bam} > {output.bamls}
        """

rule rMats:
    input:
        bam = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.rmdup.bam",
        bamls = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.bam.list",
        medRL = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.medianRL.txt",
        gtf = "resources/genome.gtf",
    output:
        se = "{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.txt",
    conda:
        "../envs/altsplice.yaml",
    wildcard_constraints:
        sample = common_constraint,
    params:
        extra = "--variable-read-length --statoff",
        # keep the directory for rmats.py to write to
        direc = "{runid}/results/rmats/{sample}/{chrom}/",
    log:
        "{runid}/logs/rmats/{sample}.{chrom}.log",
    threads: 12,
    shell:
        """
        # make sure output directory exists
        mkdir -p {params.direc}

        readLen=$(cat {input.medRL}) &&
        rmats.py --b1 {input.bamls} --readLength ${{readLen}} \
            --nthread {threads} --od {params.direc} \
            --gtf {input.gtf} --tmp {params.direc} {params.extra} \
            > {log} 2>&1

        # touch the tracked file to make Snakemake happy
        if [ ! -f {output.se} ]; then
            echo "SE.MATS.JC content placeholder" > {output.se}
        fi
        """

#rule rMats:
#    input:
#        bam = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.rmdup.bam",
#        bamls = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.bam.list",
#        #bai = "{runid}/results/reads/star/mrkdup/{sample}/{chrom}.bam.bai",
#        medRL="{runid}/results/reads/star/mrkdup/{sample}/{chrom}.medianRL.txt",
#        gtf = "resources/genome.gtf",
#    output:
#        se = "{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.txt",
#    conda:
#      "../envs/altsplice.yaml"
#    wildcard_constraints:
#        sample = common_constraint
#    params:
#        extra = "--variable-read-length --statoff",
#        direc = directory("{runid}/results/rmats/{sample}/{chrom}/"),
#        #rl = 33
#    log:
#        "{runid}/logs/rmats/{sample}.{chrom}.log",
#    threads: 12
#    shell:  
#        """
#        readLen=$(cat {input.medRL} ) &&
#        rmats.py --b1 {input.bamls} --readLength ${{readLen}} --nthread {threads} --od {params.direc} --gtf {input.gtf} --tmp {params.direc} {params.extra} > {log} 2>&1
#        """
# 


#rule grepMETse:
#    input:
#      se = "{runid}/results/rmats/{sample}/SE.MATS.JC.txt",
#      targets = config['splicingTargets']
#    output:
#      "{runid}/results/rmats/{sample}/seMET.txt",    
#    shell:
#        """
#        head -1 {input.se} > {output} &&
#        readarray -t < <(cat {input.targets})
#        for i in ${{MAPFILE[@]}}; do 
#          egrep -f {input.targets} {input.se} >> {output};
#        done
#        """

rule clinEx:
    input: 
        se = "{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.txt",
        #db = "resources/Exon_skipping_inducing_mutations_all_info.txt"
        db = "resources/exon_targets.tsv",
    output:
        clinout = "{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.clinout.tsv",
        nonclinout = "{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.nonclinout.tsv",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log: "{runid}/logs/clinEx/{sample}.{chrom}.log"
    shell:"""
        echo "exon_file = {input.se}" >&2
        python scripts/clinExSkip.py -e {input.se} -c {input.db}  2> {log}
        """


rule exonfltr:
    input:
        co = "{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.clinout.tsv",
        h5 = "resources/exons23.h5",
        bl = "resources/bl_exonSkip.tsv"
    output:
        "{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.clinout.fltrd.tsv",
    conda:
      "../envs/hts.yaml"
#    wildcard_constraints:
#        sample = common_constraint
    log: "{runid}/logs/exonfltr/{sample}/{chrom}.log"
    shell:"""
        python scripts/exonSkipfltr.py -e {input.co} -d {input.h5} -b {input.bl} > {output} 2> {log}
        """

rule cat_exonSkippers:
    input:
        expand("{runid}/results/rmats/{sample}/{chrom}/SE.MATS.JC.clinout.fltrd.tsv", runid= runid, sample= idkeys, chrom= ['chr7', 'chr17', 'chrX'])
    output:
        "{runid}/results/rmats/{sample}.exonSkipFull.tsv"
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    params:
        direc = "{runid}/results/rmats/{sample}"
    log: "{runid}/logs/cat_exonSkippers/{sample}.log"
    shell:"""
        cat {runid}/results/rmats/{wildcards.sample}/chr*/SE.MATS.JC.clinout.fltrd.tsv > {output} 2> {log}
        """


rule egfr_v3:
    input:
        bam = "{runid}/results/reads/star/mrkdup/{sample}.bam",
        bai = "{runid}/results/reads/star/mrkdup/{sample}.bam.bai",
    output:
        "{runid}/results/rmats/{sample}.egfr_v3.out"
    conda:
      "../envs/altsplice.yaml"
    log: "{runid}/logs/egfr_v3/{sample}.log"
    shell: """
        egfr-v3-determiner -r hg38 {input.bam} -w all -v all > {output} 2> {log}
        """


rule exonFinalOut:
    input:
        es = "{runid}/results/rmats/{sample}.exonSkipFull.tsv",
        egfr = "{runid}/results/rmats/{sample}.egfr_v3.out"
    output:
        "{runid}/results/rmats/{sample}.exonSkip.tsv"
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/exonFinalOut/{sample}.log"
    shell: """
        python scripts/exonFinalOut.py -e {input.es} -g {input.egfr} > {output} 2> {log}
        """



