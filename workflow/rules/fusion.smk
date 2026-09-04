
rule arriba:
    input:
        bam="{runid}/results/bam/{sample}.bam",
        genome=genome_fa,
        annotation=genome_gtf,
    output:
        fusions="{runid}/results/fusions/{sample}.fusions.tsv",
        discarded="{runid}/results/fusions/{sample}.fusions.discarded.tsv",
        #done="{runid}/results/fusions/{sample}.arriba.done",
    params:
        genome_build="GRCh38",
    log:
        "{runid}/results/fusions/{sample}.arriba.log",
    threads: 1
    shell:
        """
        DB=/opt/conda-env/var/lib/arriba
        BL=$(ls "$DB"/blacklist_hg38_GRCh38_*.tsv.gz 2>/dev/null | head -1)
        KN=$(ls "$DB"/known_fusions_hg38_GRCh38_*.tsv.gz 2>/dev/null | head -1)
        [ -n "$BL" ] && [ -n "$KN" ] || { echo "arriba DB files not found under $DB" >&2; exit 1; }
        arriba -x {input.bam} -a {input.genome} -g {input.annotation} \
            -b "$BL" -k "$KN" \
            -o {output.fusions} -O {output.discarded} -u \
            > {log} 2>&1
        """


rule arriba_draw_fusions:
    """Render fusion visualizations with arriba's draw_fusions.R.

    Draws gene structures, breakpoint positions, and (when alignments are
    provided) read-coverage panels supporting each fusion.
    """
    input:
        fusions = "{runid}/results/fusions/{sample}.fusions.tsv",
        gtf     = "resources/genome.gtf",
        bam     = "{runid}/results/bam/{sample}.bam",
    output:
        pdf = "{runid}/results/fusions/{sample}_fusions.pdf",
    log:
        "{runid}/logs/arriba_draw_fusions/{sample}.log",
    threads: 1
    shell:
        """
        draw_fusions.R \
            --fusions {input.fusions} \
            --annotation {input.gtf} \
            --output {output.pdf} \
            --alignments {input.bam} \
            > {log} 2>&1
        """


rule get_clinFuse:
    input:
        fus = "{runid}/results/arriba/{sample}/fusions.tsv",
        clinTab = "resources/winters_and_cegat_genes.tsv",
    output:
        clinout = "{runid}/results/arriba/{sample}/fusions.clinout.tsv",
        nonclinout = "{runid}/results/arriba/{sample}/fusions.nonclinout.tsv",
    conda:
      "../envs/hts.yaml"
    log: 
        "{runid}/logs/clinFuse/{sample}.log"
    shell:"""
        python scripts/clinFuse.py -f {input.fus} -c {input.clinTab} 2> {log}
        """


rule filter_arriba:
    input:
        fus = "{runid}/results/arriba/{sample}/fusions.{clin_nonclin}.tsv",
        bl = "resources/bl_exonSkip.tsv",
    output:
        "{runid}/results/arriba/{sample}/fusions.{clin_nonclin}.fltrd.tsv",
    conda:
      "../envs/hts.yaml"
    #wildcard_constraints:
    #    sample = common_constraint
    log:
        "{runid}/logs/fltr_arriba/{sample}/{clin_nonclin}.log"
    shell:"""
        python scripts/fusionfltr.py -i {input.fus} -b {input.bl} > {output} 2> {log}
        """

rule gtf2exon_h5:
    input:
      gtf = "resources/genome.gtf"
    output:
      h5 = "resources/exons_from_gtf.h5"
    conda:
      "../envs/hts.yaml"
    log:
      "logs/gtf2exon_h5.log"
    shell: """
      python scripts/gtf2exon_h5.py -g {input.gtf} -h5 {output.h5}
      """


rule getExons_arriba:
    input:
        tsv = "{runid}/results/arriba/{sample}/fusions.clinout.fltrd.tsv",
        h5 = "resources/exons_from_gtf.h5",
    output:
        "{runid}/results/arriba/{sample}/fusions.clinout.fltrd.ex.tsv",
    conda:
      "../envs/hts.yaml"
    #wildcard_constraints:
    #    sample = common_constraint
    log:
        "{runid}/logs/getExons_arriba/{sample}.log"
    shell:"""
        python scripts/fus_get_exonNo.py -i {input.tsv} -d {input.h5} > {output} 2> {log}
        """



