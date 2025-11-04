ruleorder: call_consensus_reads > fq2ubam


rule fq2ubam:
    """Generates a uBam from R1 and R2 fastq files."""
    input:
        fq1 = "{runid}/results/reads/cat/{sample}_R1.fq.gz",
        fq2 = "{runid}/results/reads/cat/{sample}_R2.fq.gz",
    output:
        bam = temp("{runid}/results/reads/{sample}.unmapped.bam")
    wildcard_constraints:
        sample= r"[^\.]+(?<!\.cons)"
    conda:
      "../envs/umi.yaml"
    resources:
        mem_gb = 14
    log:
        "{runid}/logs/fgbio/fastq_to_ubam/{sample}.log"
    shell:
        " fgbio -Xmx14g --compression 1 --async-io FastqToBam "
        "   --input {input.fq1} {input.fq2} "
        "   --read-structures +T +T "
        "   --sample {wildcards.sample} "
        "   --library {wildcards.sample} "
        "   --platform-unit flowcell.lane "
        "   --output {output.bam} &> {log} "


rule extract_umis:
    """Extract UMIs from reads into BAM tags (ZA, ZB, RX)."""
    input:
        bam = "{runid}/results/reads/{sample}.unmapped.bam"
    output:
        bam = temp("{runid}/results/reads/{sample}.umi_extracted.unmapped.bam")
    conda:
        "../envs/umi.yaml"
    wildcard_constraints:
        sample= r"[^\.]+(?<!\.cons)"
    resources:
        mem_gb = 8
    log:
        "{runid}/logs/fgbio/extract_umis/{sample}.log"
    shell:
        "fgbio -Xmx8g ExtractUmisFromBam "
        "  --input {input.bam} "
        "  --output {output.bam} "
        "  --read-structure 5M2S+T 5M2S+T "
        "  --molecular-index-tags ZA ZB "
        "  --single-tag RX "
        "  &> {log}"


rule align_bam:
    """Takes an unmapped BAM and generates an aligned BAM using bwa and ZipperBams. NOTE: this runs twice per sample, with {prefix}: 1) {sample}.umi_extracted & 2) {sample}.cons"""
    input:
        bam = "{runid}/results/reads/{prefix}.unmapped.bam",
        ref = "resources/genome.fa"
    output:
        bam = temp("{runid}/results/reads/{prefix}.mapped.bam")
    conda:
      "../envs/umi.yaml"
    wildcard_constraints:
        prefix= r"" 
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
        bam = "{runid}/results/reads/{sample}.umi_extracted.mapped.bam",
    output:
        bam = temp("{runid}/results/reads/{sample}.grouped.bam"),
        stats = "{runid}/results/reads/{sample}.grouped-family-sizes.txt"
    conda:
      "../envs/umi.yaml"
    wildcard_constraints:
        sample= r"[^\.]+(?<!\.cons)"
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
        "  --strategy adjacency "     #paired " 
        "  --edits {params.allowed_edits} "
        "  --output {output.bam} "
        "  --family-size-histogram {output.stats} "
        "  &> {log} "


rule call_consensus_reads:
    """Call consensus reads from the grouped reads using duplex strategy."""
    input:
        bam = "{runid}/results/reads/{sample}.grouped.bam",
    output:
        bam = temp("{runid}/results/reads/{sample}.cons.unmapped.bam"),
    conda:
      "../envs/umi.yaml"
    wildcard_constraints:
        sample= r"[^\.]+(?<!\.cons)"
    params:
        min_reads = 1,  # 2 1 1     # if using paired & CallDuplexConsensus
        min_base_qual = 20
    threads:
        4
    resources:
        mem_gb = 8
    log:
        "{runid}/logs/fgbio/call_consensus_reads/{sample}.log"
    shell:
        "fgbio -Xmx8g --compression 1 CallMolecularConsensusReads "     #CallDuplexConsensusReads " 
        "  --input {input.bam} "
        "  --output {output.bam} "
        "  --min-reads 1 "
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
        min_reads = 1,
        min_base_qual = 30, # 30 
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


rule bam2fq:
    input:
        #"{runid}/results/reads/{sample}.cons.mapped.bam"
        "{runid}/results/reads/{sample}.cons.filtered.bam",
    output:
        fastq1 = "{runid}/results/reads/cons/{sample}.cons.1.fq",
        fastq2 = "{runid}/results/reads/cons/{sample}.cons.2.fq",       
    conda:
      "../envs/hts.yaml"
    wildcard_constraints:
        sample = r"[^\.]+(?<!\.cons)"
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



