ruleorder: call_consensus_reads > fastq_to_ubam



rule fastq_to_ubam:
    """Generates a uBam from R1 and R2 fastq files."""
    input:
        fq1 = "{runid}/results/bcl2fq/cat/{sample}_R1.fq.gz",
        fq2 = "{runid}/results/bcl2fq/cat/{sample}_R2.fq.gz",
    params:
        rs1 = r1_read_structure,
        rs2 = r2_read_structure,
    output:
        bam = temp("{runid}/results/reads/{sample}.unmapped.bam")
    resources:
        mem_gb = 4
    log:
        "{runid}/logs/fgbio/fastq_to_ubam/{sample}.log"
    shell:
        " fgbio -Xmx4g --compression 1 --async-io FastqToBam "
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
        bam = "{runid}/results/reads/{prefix}.mapped.bam"
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
        "   | fgbio -Xmx4g --compression 1 --async-io ZipperBams "
        "       --unmapped {input.bam} "
        "       --ref {input.ref} "
        "       --output {output.bam} "
        "       --tags-to-reverse Consensus "
        "       --tags-to-revcomp Consensus "
        " ) &> {log}"
        

rule group_reads:
    """Group the raw reads by UMI and position ready for consensus calling."""
    input:
        bam = "{runid}/results/reads/{sample}.mapped.bam",
    output:
        bam = "{runid}/results/reads/{sample}.grouped.bam",
        stats = "{runid}/results/reads/{sample}.grouped-family-sizes.txt"
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
        bed = "resources/twist_rna_exome_target_regions_hg38_annotated.bed",
        # target_file_UMI_demo_data_hg38.bed",
        ref = "resources/genome.fa",
        probes = "resources/probe_file_UMI_demo_data_hg38.interval_list",
        targets = "resources/target_file_UMI_demo_data_hg38.interval_list",
    output:
        metrics = "{runid}/results/picard/metrics/cons/{sample}.metrics.tsv",
        perTargetCov = "{runid}/results/picard/metrics/cons/{sample}.cov.bed",
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
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/collectHs/cons/{sample}.formatHs.log",
    shell:
        """
        scripts/formatHs.sh {input} {runid} > {output} 2> {log}
        """

rule sambamba_markdup:
    input:
        "{runid}/results/reads/{sample}.cons.mapped.bam"
    output:
        "{runid}/results/reads/{sample}.cons.mapped.mrkdup.bam",
    priority: 20
    params:
        extra="-r"  # optional parameters
    log: "{runid}/logs/sambamba-markdup/{sample}.log"
    threads: 8
    wrapper:
        "v1.31.1/bio/sambamba/markdup"

rule bam2fq:
    input:
        "{runid}/results/reads/{sample}.cons.mapped.mrkdup.bam",
    output:
        fastq1 = "{runid}/results/reads/cons/{sample}.cons.1.fq",
        fastq2 = "{runid}/results/reads/cons/{sample}.cons.2.fq",       
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
        aln = "{runid}/results/reads/star/{sample}.bam",
        sj = "{runid}/results/reads/star/{sample}/SJ.out.tab",
    wildcard_constraints:
        sample = common_constraint
    params:
        smpl = "{sample}",
        #rg = "{rgid}",
        lib = "Library1",
        pu = "Unit1",
        pl = "Illumina",
        extra="--chimSegmentMin 12 --chimOutType WithinBAM --readFilesSAMattrKeep All --quantMode GeneCounts --limitBAMsortRAM 35000000000 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 1 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30", # --sjdbGTFfile genome.gtf ",# --outSAMattrRGline ID:{sample} --sjdbGTFfile {} {}".format(
        #"resources/genome.gtf", config["params"]["star"]
        #),
    log: "{runid}/logs/star/{sample}.log"
    threads: 24
    shell:"""
        STAR --runThreadN {threads} --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} {params.extra} --outFileNamePrefix {runid}/results/reads/star/{wildcards.sample}/ --outStd BAM_SortedByCoordinate --outSAMattrRGline ID:{rgid} SM:{params.smpl} LB:{params.pl} PU:{params.pu} > {output.aln} 
        """

rule star_markdup:
    input:
        "{runid}/results/reads/star/{sample}.bam"
    output:
        "{runid}/results/reads/star/mrkdup/{sample}.bam",
    wildcard_constraints:
        sample = common_constraint
    priority: 20
    params:
        extra=""  #"-r"  # optional parameters
    log: "{runid}/logs/star_markdup/{sample}.log"
    threads: 8
    wrapper:
        "v2.2.1/bio/sambamba/markdup"


rule filter_bam:
    input:
        "{runid}/results/reads/star/{sample}.bam",
    output:
        "{runid}/results/reads/star_fltrd/{sample}.bam",
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/fltr_bam/{sample}.log"
    shell:"""
        samtools view -h {input} | awk 'length($10) > 40 || $1 ~ /^@/' | samtools view -bS - > {output}
        """


rule arriba:
    input:
        bam="{runid}/results/reads/star_fltrd/{sample}.bam",
        genome="resources/genome.fa",
        annotation="resources/genome.gtf",
        # optional: # A custom tsv containing identified artifacts, such as read-through fusions of neighbouring genes.
        #blacklist="resources/blacklist_hg38_GRCh38_v2.3.0.tsv.gz",
        #known_fusions="resources/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz",
        # default blacklists are selected via blacklist parameter
        # see https://arriba.readthedocs.io/en/latest/input-files/#blacklist
        #custom_blacklist=[],
    output:
        fusions="{runid}/results/arriba/{sample}/fusions.tsv",
        discarded="{runid}/results/arriba/{sample}/fusions.discarded.tsv",
        #done="{runid}/results/arriba/{sample}/arriba.done",
    params:
        # required if blacklist or known_fusions is set
        genome_build="GRCh38",
        default_blacklist=True,
        default_known_fusions=True,
        extra=""  #"-u",   #alignIntronMax",
    log:
        "{runid}/logs/arriba/{sample}.log",
    threads: 1
    wrapper:
        "v1.23.4/bio/arriba"

rule filter_arriba:
    input:
        "{runid}/results/arriba/{sample}/fusions.tsv",
    output:
        "{runid}/results/arriba/{sample}/fusions.fltrd.tsv",
    #wildcard_constraints:
    #    sample = common_constraint
    log:
        "{runid}/logs/fltr_arriba/{sample}/log"
    shell:"""
        scripts/fusionfltr.py -i {input} > {output} 2> {log}
        """


rule collectHs_star:
    input:
        bam = "{runid}/results/reads/star_fltrd/{sample}.bam",
        bed = "resources/twist_rna_exome_target_regions_hg38_annotated.bed",
        # target_file_UMI_demo_data_hg38.bed",
        ref = "resources/genome.fa",
        probes = "resources/probe_file_UMI_demo_data_hg38.interval_list",
        targets = "resources/target_file_UMI_demo_data_hg38.interval_list",
    output:
        metrics = "{runid}/results/picard/metrics/star/{sample}.metrics.tsv",
        perTargetCov = "{runid}/results/picard/metrics/star/{sample}.cov.bed",
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
    wildcard_constraints:
        sample = common_constraint
    log:
        "{runid}/logs/collectHs/star/{sample}.formatHs.log",
    shell:
        """
        scripts/formatHs.sh {input} {runid} > {output} 2> {log}
        """


