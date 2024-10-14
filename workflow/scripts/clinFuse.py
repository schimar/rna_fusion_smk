
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

# This script takes as input an arriba fusions file and parses each line and respective fusion partner.If at least one of the fusion partners occurs in the list of clinically relevant fusions (in resources/winters_and_cegat_genes.tsv), then it will be output into fusions.clin.tsv. Otherwise, fusions.nonclin.tsv.

# Usage: scripts/clinFuse.py -f fusions.tsv -c ~/smk/rna_fusion_quant/workflow/resources/winters_and_cegat_genes.tsv

# -----------------------------------------------------------------------------


#lshdr = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
#hdrline = '\t'.join(lshdr)


parser = ArgumentParser()
parser.add_argument("-f", "--fusions", dest="infus", help="input fusion file name", metavar="<input>")
parser.add_argument("-c", "--clinfus", dest="clinfus", help="list of known fusion genes", metavar="<clinfus>")


#parser.add_argument("-q", "--quiet",
#                    action="store_false", dest="verbose", default=True,
#                    help="don't print status messages to stdout")

args = parser.parse_args()

fusion_file = args.infus
clinfus = args.clinfus

if not fusion_file: #.exists():
    print("Please specify a valid fusion file")
    raise SystemExit(1)

if not clinfus: #.exists():
    print("Please specify a valid file with known fusion genes")
    raise SystemExit(1)


# -------------------------------------------------------------

if __name__ == "__main__":

    with open(clinfus, 'rt') as clinf:
        # create dict from clinfus file (w/ key = gene & value = "<description>\t<fusion partner>")
        clindict = dict()
        for line in clinf:
            line = line.strip('\n')
            if line[0:2] == "G":
                continue
            else:
                lspl = line.split('\t')
                if not lspl[0] in clindict:
                    clindict[lspl[0]] = "\t".join(lspl[1:3])
                else:
                    continue

#    for key, value in clindict.items():
#        print(key, value.split('\t'))


    with open(fusion_file, 'rt') as fus:
        #split fusion_file into fusions.clinout and fusions.nonclinout, based on clindict
        ffstm = fusion_file.split('.')[0]
        clinout = open(ffstm + ".clinout.tsv", 'w')
        nonclinout = open(ffstm + ".nonclinout.tsv", 'w')
        for line in fus:
            line = line.strip('\n')
            if line[0:2] == '#g':
                clinout.write(line + f"\tgene\tname\tcommon partner\n")
                nonclinout.write(line + f"\tgene\tname\tcommon partner\n")
            elif line in ['\n', '\r\n']:
                continue
            else:
                #lspl = line.split()
                gene1, gene2 = line.split()[0:2]    #, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, confidence, reading_frame, tags, retained_protein_domains, closest_genomic_breakpoint1, closest_genomic_breakpoint2, gene_id1, gene_id2, transcript_id1, transcript_id2, direction1, direction2, filters, fusion_transcript, peptide_sequence, read_id = lspl
                #reads = list(map(int, [split_reads1, split_reads2, disco_mates, cov1, cov2]))
                if gene1 in clindict and gene2 in clindict:
                    descr = clindict[gene1].split('\t')[0], clindict[gene2].split('\t')[0]
                    newline =  line + f"\t{gene1}::{gene2}\t{descr[0]}::{descr[1]}"
                    clinout.write(newline + '\n')
                elif gene1 in clindict:
                    newline = line + f"\t{gene1}\t{clindict[gene1]}"
                    clinout.write(newline + '\n')
                elif gene2 in clindict:
                    newline = line + f"\t{gene2}\t{clindict[gene2]}"
                    clinout.write(newline + '\n')
                else:
                    nonclinout.write(line + '\n')

    fus.close()
    clinout.close()
    nonclinout.close()



