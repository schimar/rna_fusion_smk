
import glob


def _arriba_db_file(kind):
    """Return the first GRCh38 arriba database file of the given kind from the central image."""
    hits = sorted(
        glob.glob(f"/opt/conda-env/var/lib/arriba/{kind}_hg38_GRCh38_*.tsv.gz")
    )
    if not hits:
        raise FileNotFoundError(
            f"arriba {kind} database file for GRCh38 not found in /opt/conda-env/var/lib/arriba"
        )
    return hits[0]


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
        blacklist=lambda wildcards: _arriba_db_file("blacklist"),
        known_fusions=lambda wildcards: _arriba_db_file("known_fusions"),
    log:
        "{runid}/results/fusions/{sample}.arriba.log",
    threads: 1
    shell:
        """
        arriba -x {input.bam} -a {input.genome} -g {input.annotation} \
            -b {params.blacklist} -k {params.known_fusions} \
            -o {output.fusions} -O {output.discarded} -u \
            > {log} 2>&1
        """


rule arriba_draw_fusions:
    """Render fusion visualizations with arriba's draw_fusions.R.

    Draws gene structures, breakpoint positions, and supporting read counts
    per fusion. Read-coverage panels (--alignments) are disabled until the
    image ships R/GenomicAlignments (draw_fusions.R hard-fails without it);
    ideogram/circos/domains panels additionally need GenomicRanges/circlize.
    """
    input:
        fusions = "{runid}/results/fusions/{sample}.fusions.tsv",
        gtf     = genome_gtf,
    output:
        pdf = "{runid}/results/fusions/{sample}_fusions.pdf",
    log:
        "{runid}/logs/arriba_draw_fusions/{sample}.log",
    threads: 1
    shell:
        """
        draw_fusions.R \
            --fusions={input.fusions} \
            --annotation={input.gtf} \
            --output={output.pdf} \
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



