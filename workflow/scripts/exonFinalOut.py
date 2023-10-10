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
    #parser.add_argument("-q", "--quiet",
    #                    action="store_false", dest="verbose", default=True,
    #                    help="don't print status messages to stdout")

    args = parser.parse_args()
    exon_file = args.inex

    if not exon_file: #.exists():
        print("Please specify a valid exon skipping file")
        raise SystemExit(1)

    def getExonNo(enst_ense_no):
        try:
            enst_e_no = enst_ense_no.split(' ')[1].split(',')
        except:
            enst_e_no = enst_ense_no.split(',')
        exnols = list()
        enstls = list()
        ensels = list()
        for enstin in enst_e_no:
            enstinspl = enstin.split('_')
            if len(enstinspl) == 3:
                enst, ense, exno = enstinspl
                if ense == '0':
                    exnols.append(np.nan)
                    enstls.append(enst)
                    ensels.append(np.nan)
                else:
                    exnols.append(exno)
                    enstls.append(enst)
                    ensels.append(ense)
            else:
                exnols.append(np.nan)
                enstls.append(enstinspl[0])
                ensels.append(np.nan)
        exnoSer = pd.Series(exnols)
        enstSer = pd.Series(enstls)
        enseSer = pd.Series(ensels)
        match1idx = exnoSer.notna()
        exno_1st_match = exnoSer[match1idx].tolist()[0]
        enst_1st_match = enstSer[match1idx].tolist()[0]
        ense_1st_match = enseSer[match1idx].tolist()[0]
        return exno_1st_match, enst_1st_match, ense_1st_match


    #hdr = "ID", "GeneID", "geneSymbol", "chrom", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "ID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen", "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference", "TCGAabbr", "cosmic_sample"

    # ID, GeneID, geneSymbol, chrom, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, ID, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, SkipFormLen, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference, TCGAabbr, cosmic_sample
# 25

    with open(exon_file, 'rt') as exskip:
        # print new header line
        header = '\t'.join(["rmatsID", "geneID", "geneSymbol", "chrom", "strand", "exonStart_0base", "exonEnd", "exonNo", "IJC", "SJC", "IncLevel", "tcga_code", "tcga_sample", "enstMatch", "enseMatch"])
        print(header)
        for line in exskip:
            line = line.strip('\n')
            lspl = line.split('\t')
            #print(lspl)
            ID, geneID, geneSymbol, chrom, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, ID, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, SkipFormLen, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference, tcga_code, tcga_sample, g1, exonSkipDB_ensts, chrom2, exPos, d1, d2, d3, d4, site, enst_ense_no = lspl[0:35]
            geneSymbol = geneSymbol.split('"')[1]
            exonNo, enstMatch, enseMatch = getExonNo(enst_ense_no)
            outline = [ID, geneID, geneSymbol, chrom, strand, exonStart_0base, exonEnd, exonNo, IJC_SAMPLE_1, SJC_SAMPLE_1, IncLevel1, tcga_code, tcga_sample, enstMatch, enseMatch]
            print('\t'.join(outline))


