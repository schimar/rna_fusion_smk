ruleorder: call_consensus_reads > fq2ubam



rule fq2ubam:
    """Generates a uBam from R1 and R2 fastq files."""
    input:
        fq1 = "{runid}/results/bcl2fq/cat/{sample}_R1.fq.gz",
        fq2 = "{runid}/results/bcl2fq/cat/{sample}_R2.fq.gz",
        #fq1 = "{runid}/results/bcl2fq/{sample}_R1_001.fastq.gz",
        #fq2 = "{runid}/results/bcl2fq/{sample}_R2_001.fastq.gz",

        #fq1 = "/mnt/sda/rnaSeq/runs/231025/results/bcl2fq/cat/{sample}_R1.fq.gz",
        #fq2 = "/mnt/sda/rnaSeq/runs/231025/results/bcl2fq/cat/{sample}_R2.fq.gz",
    params:
        rs1 = r1_read_structure,
        rs2 = r2_read_structure,
    output:
        bam = temp("{runid}/results/reads/{sample}.unmapped.bam")
    conda:
      "../envs/umi.yaml"
    resources:
        mem_gb = 14
    log:
        "{runid}/logs/fgbio/fastq_to_ubam/{sample}.log"
    shell:
        " fgbio -Xmx14g --compression 1 --async-io FastqToBam "
        "   --input {input.fq1} {input.fq2} "
        "   --read-structures {params.rs1} {params.rs2} "
        "   --sample {wildcards.sample} "
        "   --library {wildcards.sample} "
        "   --platform-unit flowcell.lane "
        "   --output {output.bam} &> {log} "


rule align_bam:
    """Takes an unmapped BAM and generates an aligned BAM using bwa and ZipperBams."""
    input:
        bam = "{runid}/results/reads/{prefix}.unmapped.bam",
        ref = "resources/genome.fa"
    output:
        bam = temp("{runid}/results/reads/{prefix}.mapped.bam")
    conda:
      "../envs/umi.yaml"
    threads:
        16
    resources:
        mem_gb = 14
    log:
        "{runid}/logs/bwa/align_bam/{prefix}.log"
    shell:
        " ( "
        " samtools fastq {input.bam} "
        "   | bwa mem -t {threads} -p -K 150000000 -Y {input.ref} - "
        "   | fgbio -Xmx8g --compression 1 --async-io ZipperBams "
        "       --unmapped {input.bam} "
        "       --ref {input.ref} "
        "       --output {output.bam} "
        "       --tags-to-reverse Consensus "
        "       --tags-to-revcomp Consensus "
        " ) &> {log}"
        
          
rule flagstat_mapped:
    input:
        "{runid}/results/reads/{sample}.mapped.bam",
    output:
        "{runid}/results/reads/fstat/{sample}.mapped.fstat",
    conda:
      "../envs/umi.yaml"
    shell:"""
        samtools flagstat {input} > {output}
        """
      

rule group_reads:
    """Group the raw reads by UMI and position ready for consensus calling."""
    input:
        bam = "{runid}/results/reads/{sample}.mapped.bam",
    output:
        bam = temp("{runid}/results/reads/{sample}.grouped.bam"),
        stats = "{runid}/results/reads/{sample}.grouped-family-sizes.txt"
    conda:
      "../envs/umi.yaml"
    params:
        allowed_edits = 1,
    threads:
        2
    resources:
        mem_gb = 12
    log:
        "{runid}/logs/fgbio/group_reads/{sample}.log"
    shell:
        "fgbio -Xmx12g --compression 1 --async-io GroupReadsByUmi "
        "  --input {input.bam} "
        "  --strategy Adjacency "
        "  --edits {params.allowed_edits} "
        "  --output {output.bam} "
        "  --family-size-histogram {output.stats} "
        "  &> {log} "


rule call_consensus_reads:
    """Call consensus reads from the grouped reads."""
    input:
        bam = "{runid}/results/reads/{sample}.grouped.bam",
    output:
        bam = temp("{runid}/results/reads/{sample}.cons.unmapped.bam"),
    conda:
      "../envs/umi.yaml"
    params:
        min_reads = 1,
        min_base_qual = 20
    threads:
        4
    resources:
        mem_gb = 8
    log:
        "{runid}/logs/fgbio/call_consensus_reads/{sample}.log"
    shell:
        "fgbio -Xmx8g --compression 1 CallMolecularConsensusReads "
        "  --input {input.bam} "
        "  --output {output.bam} "
        "  --min-reads {params.min_reads} "
        "  --min-input-base-quality {params.min_base_qual} "
        "  --threads {threads} "
        "  &> {log}"


rule filter_consensus_reads:
    """Filters the consensus reads and then sorts into coordinate order."""
    input:
        bam = "{runid}/results/reads/{sample}.cons.mapped.bam",
        ref = "resources/genome.fa",
    output:
        bam = "{runid}/results/reads/{sample}.cons.filtered.bam",
    conda:
      "../envs/umi.yaml"
    params:
        min_reads = 3,
        min_base_qual = 40,
        max_error_rate = 0.2
    threads:
        8
    resources:
        mem_gb = 8
    log:
        "{runid}/logs/fgbio/filter_consensus_reads.{sample}.log"
    shell:
        " ( "
        " fgbio -Xmx8g --compression 0 FilterConsensusReads "
        "   --input {input.bam} "
        "   --output /dev/stdout "
        "   --ref {input.ref} "
        "   --min-reads {params.min_reads} "
        "   --min-base-quality {params.min_base_qual} "
        "   --max-base-error-rate {params.max_error_rate} "
        "   | samtools sort --threads {threads} -o {output.bam}##idx##{output.bam}.bai --write-index "
        " ) &> {log} "


rule collectHs_cons:
    input:
        bam = "{runid}/results/reads/{sample}.cons.mapped.bam",
        bed = config['bed'],		#"resources/twist_rna_exome_target_regions_hg38_annotated.bed",
        # target_file_UMI_demo_data_hg38.bed",
        ref = "resources/genome.fa",
        probes = "resources/probe_file_UMI_demo_data_hg38.interval_list",
        targets = "resources/target_file_UMI_demo_data_hg38.interval_list",
    output:
        metrics = "{runid}/results/picard/metrics/cons/{sample}.metrics.tsv",
        perTargetCov = "{runid}/results/picard/metrics/cons/{sample}.cov.bed",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/collectHs/cons/{sample}.log",
    resources:
        mem_gb=4,
    threads: 8
    shell:
        """
        picard CollectHsMetrics -Xmx4g -I {input.bam} -O {output.metrics} -R {input.ref} --BAIT_INTERVALS {input.probes} --TARGET_INTERVALS {input.probes} --PER_TARGET_COVERAGE {output.perTargetCov} 2> {log}
        """


rule formatHs_cons:
    input:
        "{runid}/results/picard/metrics/cons/{sample}.metrics.tsv",
    output:
        "{runid}/results/picard/metrics/cons/{sample}.metricsFrmt.tsv",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/collectHs/cons/{sample}.formatHs.log",
    shell:
        """
        scripts/formatHs.sh {input} {runid} > {output} 2> {log}
        """

# for now, we'll omit removing duplicates at this stage, since it'll be done after star alignment
rule sambamba_markdup:
    input:
        "{runid}/results/reads/{sample}.cons.mapped.bam"
    output:
        "{runid}/results/reads/{sample}.cons.mapped.mrkdup.bam",
    conda:
      "../envs/hts.yaml"
    priority: 20
    params:
        extra="-r"  # optional parameters
    log: "{runid}/logs/sambamba-markdup/{sample}.log"
    threads: 8
    wrapper:
        "v1.31.1/bio/sambamba/markdup"

rule bam2fq:
    input:
        "{runid}/results/reads/{sample}.cons.mapped.bam",
    output:
        fastq1 = temp("{runid}/results/reads/cons/{sample}.cons.1.fq"),
        fastq2 = temp("{runid}/results/reads/cons/{sample}.cons.2.fq"),
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint,
    log:
        "{runid}/logs/samtools/bam2fq/{sample}.log",
    params:
        sort = "-m 4G",
        fastq="-n",
    threads: 8
    wrapper:
        "v2.1.1/bio/samtools/fastq/separate"



rule map_star:
    input:
        fq1 = "{runid}/results/reads/cons/{sample}.cons.1.fq",
        fq2 = "{runid}/results/reads/cons/{sample}.cons.2.fq",       
        idx="resources/star_genome",
        gtf = "resources/genome.gtf",
    output:
        aln = temp("{runid}/results/reads/star/{sample}.bam"),
        sj = "{runid}/results/reads/star/{sample}/SJ.out.tab",
    conda:
      "../envs/star.yaml"
    wildcard_constraints:
        sample = common_constraint
    params:
        smpl = "{sample}",
        #rg = "{rgid}",
        lib = "Library1",
        pu = "Unit1",
        pl = "Illumina",
        extra="--chimSegmentMin 12 --chimOutType WithinBAM --readFilesSAMattrKeep All --quantMode GeneCounts --limitBAMsortRAM 110000000000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30,  --sjdbGTFfile resources/genome.gtf"  # --outSAMattrRGline ID:{sample} --sjdbGTFfile {} {}".format(
        #"resources/genome.gtf", config["params"]["star"]
        #),
    log: "{runid}/logs/star/{sample}.log"
    threads: 24
    shell:"""
        STAR --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} {params.extra} --outFileNamePrefix {runid}/results/reads/star/{wildcards.sample}/ --outStd BAM_SortedByCoordinate --outSAMattrRGline ID:{rgid} SM:{params.smpl} LB:{params.pl} PU:{params.pu} > {output.aln} 2> {log}
        """

rule star_index_dup:
    input:
        "{runid}/results/reads/star/{sample}.bam",
    output:
        "{runid}/results/reads/star/{sample}.bam.bai",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log: "{runid}/logs/star_index_dup/{sample}.log"
    shell: """
        samtools index {input}
        """


rule star_markdup:
    input:
        "{runid}/results/reads/star/{sample}.bam"
    output:
        "{runid}/results/reads/star/mrkdup/{sample}.bam",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    priority: 20
    params:
        extra=""  #"-r"  # optional parameters
    log: "{runid}/logs/star_markdup/{sample}.log"
    threads: 8
    wrapper:
        "v2.2.1/bio/sambamba/markdup"

rule picard_markdup:
    input:
      "{runid}/results/reads/star/{sample}.bam"
    output:
      bam = temp("{runid}/results/reads/star/mrkdup/picard/{sample}.bam"),
      metrics = "{runid}/results/reads/star/mrkdup/picard/{sample}.txt"
    conda:
      "../envs/hts.yaml"
    params:
      rn_re = "'^([A-Z0-9]+):([0-9]+):([A-Z0-9]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) .*'"
    log: "{runid}/logs/picard/markdup/{sample}.log"
    threads: 2
    shell:"""
      picard MarkDuplicates I={input} O={output.bam} M={output.metrics} READ_NAME_REGEX={params.rn_re} OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 2> {log}
      """


#rule filter_bam:
#    input:
#        "{runid}/results/reads/star/{sample}.bam",
#    output:
#        "{runid}/results/reads/star_fltrd/{sample}.bam",
#    wildcard_constraints:
#        sample = common_constraint
#    log:
#        "{runid}/logs/fltr_bam/{sample}.log"
#    shell:"""
#        samtools view -h {input} | awk 'length($10) > 40 || $1 ~ /^@/' | samtools view -bS - > {output}
#        """

rule tpmCalc:
    input:
        bam = "{runid}/results/reads/star/mrkdup/{sample}.bam",
        gtf = "resources/genome.gtf",
    output:
        "{runid}/results/reads/star/mrkdup/tpm/{sample}_genes.out",
    conda:
      "../envs/hts.yaml"
    params:
	    outdir = "{runid}/results/reads/star/mrkdup/tpm"
    log:
        "{runid}/logs/TPMcalc/{sample}.log"
    shell:"""
        TPMCalculator -g {input.gtf} -b {input.bam} -o {params.outdir} 2> {log}
        """
 

rule collectHs_star:
    input:
        bam = "{runid}/results/reads/star/{sample}.bam",
        bed = config['bed'],		#"resources/twist_rna_exome_target_regions_hg38_annotated.bed",
        # target_file_UMI_demo_data_hg38.bed",
        ref = "resources/genome.fa",
        probes = "resources/probe_file_UMI_demo_data_hg38.interval_list",
        targets = "resources/target_file_UMI_demo_data_hg38.interval_list",
    output:
        metrics = "{runid}/results/picard/metrics/star/{sample}.metrics.tsv",
        perTargetCov = "{runid}/results/picard/metrics/star/{sample}.cov.bed",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/collectHs/star/{sample}.log",
    resources:
        mem_gb=8,
    threads: 8
    shell:
        """
        picard CollectHsMetrics -Xmx4g -I {input.bam} -O {output.metrics} -R {input.ref} --BAIT_INTERVALS {input.probes} --TARGET_INTERVALS {input.probes} --PER_TARGET_COVERAGE {output.perTargetCov} 2> {log}
        """


rule formatHs_star:
    input:
        "{runid}/results/picard/metrics/star/{sample}.metrics.tsv",
    output:
        "{runid}/results/picard/metrics/star/{sample}.metricsFrmt.tsv",
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/collectHs/star/{sample}.formatHs.log",
    shell:
        """
        scripts/formatHs.sh {input} {runid} > {output} 2> {log}
        """



