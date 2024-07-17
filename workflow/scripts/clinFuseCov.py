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

def gini_coefficient(coverage_values):
    sorted_values = np.sort(coverage_values)
    n = len(coverage_values)
    cumulative_values = np.cumsum(sorted_values)
    return (2.0 / n) * (sum((i + 1) * cumulative_values[i] for i in range(n))) / cumulative_values[-1] - (n + 1) / n


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

    #print('\t'.join(['gene', 'n10', 'mean', 'median', 'n90', 'prop10x', 'gini_coeff']))
    print('\t'.join(['gene', 'n10', 'prop_bases_over_10x']))
    for gene, cov in genedict.items():
        cov = [int(i) for i in cov]
        prop10x = np.round(sum(np.array(cov) > 10) / len(cov), 2)
        #n90 = int(np.percentile(cov, 90))
        n10 = int(np.percentile(cov, 10))
        #gini = np.round(gini_coefficient(cov), 2)
        #mean = np.round(np.mean(cov))
        #med = np.median(cov)
        #outls = [gene, str(n10), str(mean), str(med), str(n90), str(prop10x*100), str(gini)]
        outls = [gene, str(n10), str(prop10x*100)]
        print('\t'.join(outls))



