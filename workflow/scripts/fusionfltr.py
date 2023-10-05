#! /usr/bin/python3


# ---------------------------

#import os
#import glob
import sys
from sys import argv
from argparse import ArgumentParser
import os.path


import numpy as np
#np.random.seed(42)
import pandas as pd

# Usage: scripts/fusionfltr.py -i fusions.tsv

# -----------------------------------------------------------------------------

# filter definitions
#spl_read_thresh = 0
#cov_thresh = 0.01


# -----------------------------------------------------------------------------



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="infus",
                        help="input fusion file name", metavar="<input>")
    #parser.add_argument("-q", "--quiet",
    #                    action="store_false", dest="verbose", default=True,
    #                    help="don't print status messages to stdout")

    args = parser.parse_args()
    fusion_file = args.infus

    if not fusion_file: #.exists():
        print("Please specify a valid fusion file")
        raise SystemExit(1)


    with open(fusion_file, 'rt') as fus:
        for line in fus:
            line = line.strip('\n')
            if line[0:2] == '#g':
                if len(line.split()) <= 33:
                    header = '\t'.join(["gene1", "gene2", "strand1", "strand2", "breakpoint1", "breakpoint2", "site1", "site2", "typ", "split_reads1", "split_reads2", "disco_mates", "cov1", "cov2", "confidence", "reading_frame", "gene_id1", "gene_id2", "transcript_id1", "transcript_id2", "filters", "s2_gene", "s2_descr"])
                else:
                    header = '\t'.join(["gene1", "gene2", "strand1", "strand2", "breakpoint1", "breakpoint2", "site1", "site2", "typ", "split_reads1", "split_reads2", "disco_mates", "cov1", "cov2", "warning", "confidence", "reading_frame", "gene_id1", "gene_id2", "transcript_id1", "transcript_id2", "filters", "s2_gene", "s2_descr"])
                print(header)
            else:
                lspl = line.split()
                reading_frame = lspl[15]
                if len(lspl) > 33:
                    warn = ""
                    gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, confidence, reading_frame, tags, retained_protein_domains, closest_genomic_breakpoint1, closest_genomic_breakpoint2, gene_id1, gene_id2, transcript_id1, transcript_id2, direction1, direction2, filters, fusion_transcript, peptide_sequence, read_id, s2_gene, s2_descr = lspl[0:32]
                else:
                    gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, confidence, reading_frame, tags, retained_protein_domains, closest_genomic_breakpoint1, closest_genomic_breakpoint2, gene_id1, gene_id2, transcript_id1, transcript_id2, direction1, direction2, filters, fusion_transcript, peptide_sequence, read_id = lspl[0:30]
                split_reads1, split_reads2, disco_mates, cov1, cov2 = list(map(int, [split_reads1, split_reads2, disco_mates, cov1, cov2]))
                cov = int(cov1+cov2)
                supp_reads = int(split_reads1 + split_reads2 + disco_mates)
                warn = ""
                if len(lspl) > 33:
                    newlist = [gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, warn, confidence, reading_frame, gene_id1, gene_id2, transcript_id1, transcript_id2, filters, s2_gene, s2_descr]
                    newline = '\t'.join(map(str, newlist))
                else:
                    newlist = [gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, warn, confidence, reading_frame, gene_id1, gene_id2, transcript_id1, transcript_id2, filters]
                    newline = '\t'.join(map(str, newlist))
                if confidence != 'low':
                    if split_reads1 == 0 or split_reads2 == 0:
                        if disco_mates == 0:
                            continue
                        else:
                            if supp_reads < 5 and supp_reads >=3:
                                warn = "low supp. reads"
                                newlist[14] = warn
                                newline = '\t'.join(map(str, newlist))
                                print(newline)
                            elif supp_reads < 3:
                                continue
                    else:
                        if supp_reads < 5 and supp_reads >=3:
                            warn = "low supp. reads"
                            newlist[14] = warn
                            newline = '\t'.join(map(str, newlist))
                            print(newline)
                        elif supp_reads < 3:
                            continue
                        else:
                            print(newline)

    fus.close()

# see https://arriba.readthedocs.io/en/latest/output-files/ for details on the output
# see https://arriba.readthedocs.io/en/latest/internal-algorithm/ for details on filters

#gene1	gene2	strand1(gene/fusion)	strand2(gene/fusion)	breakpoint1	breakpoint2	site1	site2	type	split_reads1	split_reads2	discordant_mates	coverage1	coverage2	confidence	reading_frame	tags	retained_protein_domains	closest_genomic_breakpoint1	closest_genomic_breakpoint2	gene_id1	gene_id2	transcript_id1	transcript_id2	direction1	direction2	filters	fusion_transcript	peptide_sequence	read_identifiers


# -------------------------------------------------------------------------

    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
#   chr129446335.A<INV>..TYPE=fusion;CHR2=chr1;END=42491871;PE=7359
# Example: A fusion was detected when comparing the patient sequence and the reference genome GRCh37. The first break point is estimated to be on chromosome 1 at position 29446335 and the second break point is approximately at position 42491871 on chromosome 2. This variant is supported by 7359 reads.
