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


# Usage: scripts/clinExSkip.py -e exons.tsv -c ~/smk/rna_fusion_quant/workflow/resources/exon_targets.tsv

# -----------------------------------------------------------------------------



parser = ArgumentParser()
parser.add_argument("-e", "--ex_file", dest="inex", help="input exon skipping file name", metavar="<input>")
parser.add_argument("-c", "--clinex", dest="clinex", help="list of known exon skipping genes", metavar="<clinex>")


#parser.add_argument("-q", "--quiet",
#                    action="store_false", dest="verbose", default=True,
#                    help="don't print status messages to stdout")

args = parser.parse_args()

exon_file = args.inex
clinex = args.clinex

#if not exon_file: #.exists():
#    print("Please specify a valid exon skipping file")
#    raise SystemExit(1)

if not clinex: #.exists():
    print("Please specify a valid file with known exon skipping genes")
    raise SystemExit(1)


# -------------------------------------------------------------

if __name__ == "__main__":

    with open(clinex, 'rt') as clinf:
        # create dict from clinex file (w/ key = gene & value = "<description>\t<fusion partner>")
        clindict = dict()
        for line in clinf:
            line = line.strip('\n')
            if line[0:2] == "G":
                continue
            else:
                lspl = line.split('\t')
                if not lspl[2] in clindict:
                    clindict[lspl[2]] = "\t".join(lspl[0:12])
                else:
                    continue

#    for key, value in clindict.items():
#        print(key, value.split('\t'))


    with open(exon_file, 'rt') as exskip:
        #split exon_file into fusions.clinout and fusions.nonclinout, based on clindict
        ffstm = ".".join(exon_file.split('.')[0:3])
        #print(".".join(ffstm))
        clinout = open(ffstm + ".clinout.tsv", 'w')
        nonclinout = open(ffstm + ".nonclinout.tsv", 'w')
        for line in exskip:
            line = line.strip('\n')
            if line[0:2] == 'ID':
                #print(line)
                clinout.write(line + f"\tTCGAabbr\tcosmic_sample\t\t\t\t\t\t\t\t\t\t\t\t\n")
                nonclinout.write(line + f"\tTCGAabbr\tcosmic_sample\t\t\t\t\t\t\t\t\t\t\t\t\n")
            elif line in ['\n', '\r\n']:
                continue
            else:
                lspl = line.split()
                gene = lspl[2].strip('"')
                #print(gene)
                if gene in clindict:
                    cdg = clindict[gene]
                    newline = line + f"\t{cdg}"
                    clinout.write(newline + '\n')
                else:
                    #cdg = clindict[gene]
                    #newline = line + f"{cdg}"
                    nonclinout.write(line + '\n')
                #gene1, line.split()[0:2]    #, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, confidence, reading_frame, tags, retained_protein_domains, closest_genomic_breakpoint1, closest_genomic_breakpoint2, gene_id1, gene_id2, transcript_id1, transcript_id2, direction1, direction2, filters, fusion_transcript, peptide_sequence, read_id = lspl
                #reads = list(map(int, [split_reads1, split_reads2, disco_mates, cov1, cov2]))
#                if gene1 in clindict and gene2 in clindict:
#                    descr = clindict[gene1].split('\t')[0], clindict[gene2].split('\t')[0]
#                    newline =  line + f"\t{gene1}::{gene2}\t{descr[0]}::{descr[1]}"
#                    clinout.write(newline + '\n')
#                elif gene1 in clindict:
#                    newline = line + f"\t{gene1}\t{clindict[gene1]}"
#                    clinout.write(newline + '\n')
#                elif gene2 in clindict:
#                    newline = line + f"\t{gene2}\t{clindict[gene2]}"
#                    clinout.write(newline + '\n')
#                else:
#                    nonclinout.write(line + '\n')

    exskip.close()
#    clinout.close()
#    nonclinout.close()



