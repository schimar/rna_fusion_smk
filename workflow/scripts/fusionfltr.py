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

#import warnings
#warnings.filterwarnings('ignore')


#if len(sys.argv) == 2:
#    if os.path.exists(sys.argv[1]):
#        with open(sys.argv[1]) as file:
#            print(file.read())
#    else:
#        print('No such file')
#elif len(sys.argv) < 2:
#    print('Too few arguments')
#else:
#    print('Too many arguments')

# -----------------------------------------------------------------------------

# filter definitions
#spl_read_thresh = 0
#cov_thresh = 0.01


# -----------------------------------------------------------------------------


#lshdr = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
#hdrline = '\t'.join(lshdr)


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

# NOTE: move the argparser below if name == main
# ---------------------------------- #

if __name__ == "__main__":

    with open(fusion_file, 'rt') as fus:
        for line in fus:
            line = line.strip('\n')
            if line[0:2] == '#g':
                #continue
                print(line)
            else:
                lspl = line.split()
                reading_frame = lspl[15]
                gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, confidence, reading_frame, tags, retained_protein_domains, closest_genomic_breakpoint1, closest_genomic_breakpoint2, gene_id1, gene_id2, transcript_id1, transcript_id2, direction1, direction2, filters, fusion_transcript, peptide_sequence, read_id = lspl[0:30]
                #split_reads1, split_reads2, disco_mates, cov1, cov2 = list(map(int, [split_reads1, split_reads2, disco_mates, cov1, cov2]))
                #expr_by_cov = int(split_reads1+split_reads2+disco_mates) / int(cov1+cov2)
                #if expr_by_cov >= 0.01 and confidence != "low":
               # if gene1 == 'SLC45A3' and gene2 == 'BRAF':
                    # this one doesn't show up as in-frame nor out-of-frame (5'UTR/splice-site)
               #     print(line)
                if confidence != 'low':
                    if reading_frame == 'in-frame' or reading_frame == '.':
                        print(line)
                    #continue
                #else:
                #    print(split_reads1, split_reads2, disco_mates, cov1, cov2, expr_by_cov, confidence)#expr_by_cov, confidence)


                #if (int(split_reads1) >= spl_read_thresh and int(split_reads2) >= spl_read_thresh and int(cov1) >= cov_thresh and int(cov2) >= cov_thresh):
                    # get chroms for both from breakpoints

                    #info =
                    #newline = '\t'.join(



                    #print(breakpoint1, breakpoint2, split_reads1, split_reads2, cov1, cov2, gene1, gene2, filters)


    fus.close()

# see https://arriba.readthedocs.io/en/latest/output-files/ for details on the output
# see https://arriba.readthedocs.io/en/latest/internal-algorithm/ for details on filters

#gene1	gene2	strand1(gene/fusion)	strand2(gene/fusion)	breakpoint1	breakpoint2	site1	site2	type	split_reads1	split_reads2	discordant_mates	coverage1	coverage2	confidence	reading_frame	tags	retained_protein_domains	closest_genomic_breakpoint1	closest_genomic_breakpoint2	gene_id1	gene_id2	transcript_id1	transcript_id2	direction1	direction2	filters	fusion_transcript	peptide_sequence	read_identifiers


# -------------------------------------------------------------------------

    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
#   chr129446335.A<INV>..TYPE=fusion;CHR2=chr1;END=42491871;PE=7359
# Example: A fusion was detected when comparing the patient sequence and the reference genome GRCh37. The first break point is estimated to be on chromosome 1 at position 29446335 and the second break point is approximately at position 42491871 on chromosome 2. This variant is supported by 7359 reads.
