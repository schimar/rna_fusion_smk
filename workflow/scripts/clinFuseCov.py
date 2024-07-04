#! /usr/bin/python3


# ---------------------------

#import os
#import glob
#import sys
#from sys import argv
from argparse import ArgumentParser
#import os.path
#import h5py


import numpy as np
#np.random.seed(42)
#import pandas as pd

# Usage: scripts/clinFuseCov.py -c bcov.tsv -d resources/s2_winters2018.tsv

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--cov_file", dest="incov",
                        help="input bedtools coverage file name", metavar="<input>")
    parser.add_argument("-d", "--clinData", dest="winters",
                        help="clinical targets file name", metavar="<h5>")
    args = parser.parse_args()


    if not args.incov: #.exists():
        print("Please specify a valid bedtools coverage file")
        raise SystemExit(1)
    if not args.winters:
        print("Please specify a valid clinical targets file name")
        raise SystemExit(1)



    with open(args.winters, 'rt') as winters:
        targets = list()
        for line in winters:
            if line[0:2] != '#G':
                gene = line.split('\t')[0]
                targets.append(gene)


    with open(args.incov, 'rt') as cov:
        genedict = {}
        for line in cov:
            lspl = line.strip('\n').split('\t')
            gene, pos, cov = lspl[3:6]
            if gene in targets:
                if not gene in genedict:
                    genedict[gene] = [cov]
                else:
                    genedict[gene].append(cov)


    for gene, cov in genedict.items():
        cov = [int(i) for i in cov]
        n10 = int(np.round(np.percentile(cov, 10)))
        print('\t'.join([gene, str(n10)]))



