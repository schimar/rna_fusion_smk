#! /usr/bin/python3


# ---------------------------

#import os
#import glob
import sys
from sys import argv
from argparse import ArgumentParser
import os.path
import h5py
from codecs import decode


import numpy as np
#np.random.seed(42)
import pandas as pd

# Usage: scripts/exonSkipfltr.py -e SE.MATS.JC.clinout.tsv -d resources/exons23.h5

# -----------------------------------------------------------------------------

# filter definitions
norm_ratio_th = 0.5
sjc_th = 40
strandls = ['+']

# -----------------------------------------------------------------------------



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-e", "--ex_file", dest="inex",
                        help="input exon skipping file name", metavar="<input>")
    parser.add_argument("-d", "--h5db", dest="h5in",
                        help="hdf5 file name", metavar="<h5>")
    parser.add_argument("-b", "--blacklist", dest="blackls",
                        help="input blacklist file name", metavar="<blackls>")

    #parser.add_argument("-q", "--quiet",
    #                    action="store_false", dest="verbose", default=True,
    #                    help="don't print status messages to stdout")

    args = parser.parse_args()
    exon_file = args.inex
    h5 = args.h5in
    bls = args.blackls

    if not exon_file: #.exists():
        print("Please specify a valid exon skipping file")
        raise SystemExit(1)
    if not h5:
        print("Please specify a valid hdf5 db file name")
        raise SystemExit(1)
    if not bls: #.exists():
        print("Please specify a valid blacklist file")
        raise SystemExit(1)


    def getENSTdf(f, enst):
        df = pd.DataFrame(f[enst][()])
        for col, dtype in df.dtypes.items():
            if dtype == object:  # only process bytes object columns
                df[col] = df[col].apply(lambda x: x.decode('utf-8'))
        return df


    def getExonPosLS(f, enstls, start):
        start_ens = start + 1
        xon = None
        xls = list()
        for enst in enstls:
            enst = enst.split('.')[0]
            if enst not in f:
                xls.append(f"{enst}_NA")
            else:
                xon = getENSTdf(f, enst)
                xon_start = xon[xon.exon_chrom_start == start_ens]
                if xon_start.size != 0:
                    xon_no = xon_start.index[0] + 1
                    ense = np.unique(xon[xon.exon_chrom_start == start_ens].ensembl_exon_id).astype('|S')[0].decode('utf-8')
                    xls.append(f"{enst}_{ense}_{xon_no}")
                else:
                    xls.append(f"{enst}_0_NA")
        return ','.join(xls)



    with open(bls, 'rt') as bl:
        bldict = dict()
        for line in bl:
            line = line.strip('\n')
            lspl  = line.split('\t')
            gene, chrom, strand, exonStart, exonEnd = lspl[0:5]
            if gene not in bldict:
                bldict[gene] = '_'.join([exonStart, exonEnd])

#for key, val in bldict.items():
#    print(key, val)


    with h5py.File(h5, 'r') as f:
        with open(exon_file, 'rt') as exskip:
            exon_id_ls = list()     # add id to this list, if the skipping event is accepted & printed
            for line in exskip:
                line = line.strip('\n')
                if line[0:2] == 'ID':
                    header = line
                    #print(line)
                else:
                    lspl = line.split('\t')
                    exon_id = lspl[0]
                    gene = lspl[2].strip('"')
                    strand = lspl[4]
                    ijc = int(lspl[12])
                    sjc = int(lspl[13])
                    enstls = lspl[26].split(',')
                    start = int(lspl[5])
                    stop = int(lspl[6])
                    if lspl[20] == 'NA':
                        continue
                    else:
                        if strand in strandls:
                            norm_ratio = float(lspl[20])
                            if sjc == 0:
                                continue
                            else:
                                if gene in bldict and '_'.join(map(str, [start, stop])) == bldict[gene]:
                                    continue
                                else:
                                    if norm_ratio < norm_ratio_th:
                                        if sjc >= sjc_th:
                                            exonPos = getExonPosLS(f, enstls, start)
                                            print(line, exonPos)


    #hdr = "ID      GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES       downstreamEE    ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen        PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference      TCGAabbr        cosmic_sample"


