#! /usr/bin/python3


# ---------------------------

#import os
#import glob
import sys
from sys import argv
from argparse import ArgumentParser
import os.path

import h5py
import numpy as np
#np.random.seed(42)
import pandas as pd

# Usage: scripts/fusGetExons.py -i fusions.tsv -d

# -----------------------------------------------------------------------------

# filter definitions
#spl_read_thresh = 0
#cov_thresh = 0.01


# -----------------------------------------------------------------------------



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="infus",
                        help="input fusion file name", metavar="<input>")
    parser.add_argument("-d", "--h5db", dest="h5in",
                        help="hdf5 file name", metavar="<h5>")

    #parser.add_argument("-q", "--quiet",
    #                    action="store_false", dest="verbose", default=True,
    #                    help="don't print status messages to stdout")

    args = parser.parse_args()
    fusion_file = args.infus
    h5 = args.h5in

    if not fusion_file: #.exists():
        print("Please specify a valid fusion file")
        raise SystemExit(1)

    if not h5:
        print("Please specify a valid hdf5 db file name")
        raise SystemExit(1)

    def getENSTdf(f, enst):
        #if enst not in f:
            #continue
        #    print('')#ENST not in database')
        #else:
        df = pd.DataFrame(f[enst][()])
        for col, dtype in df.dtypes.items():
            if dtype == object:  # only process bytes object columns
                df[col] = df[col].apply(lambda x: x.decode('utf-8'))
        return df

    with h5py.File(h5, 'r') as f:
        with open(fusion_file, 'rt') as fus:
            for line in fus:
                line = line.strip('\n')
                if line[0:2] == 'ge':
                    print(line)
                    header = '\t'.join(["gene1", "gene2", "strand1", "strand2", "breakpoint1", "breakpoint2", "site1", "site2", "typ", "split_reads1", "split_reads2", "disco_mates", "cov1", "cov2", "warning", "confidence", "reading_frame", "gene_id1", "gene_id2", "transcript_id1", "transcript_id2", "filters", "s2_gene", "s2_descr"])
                    #print(len(header.split('\t')))
                else:
                    lspl = line.split('\t')
                    [gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, warning, confidence, reading_frame, ensg1, ensg2, enst1, enst2, filters, s2_gene, s2_descr] = lspl
                    ensts = lspl[19:21]
                    enstls = [x.split('.')[0] for x in ensts]
                    breakpoints = [int(x.split(':')[1]) for x in [breakpoint1, breakpoint2]]
                    genels = [gene1, gene2]
                    outls = list()
                    for i, enst in enumerate(enstls):
                        try:
                            xon = getENSTdf(f, enst)
                        except:
                            xon = pd.DataFrame()
                        if xon.empty:
                            outls.append('_'.join(['NA', genels[i]]))
                        else:
                            xon_match = xon[(xon.exon_chrom_start <= breakpoints[i]) & (xon.exon_chrom_end >= breakpoints[i])]
                            if xon_match.empty:
                                outls.append('_'.join(['NA', genels[i]]))
                            else:
                                xon_no = xon_match.index[0] + 1
                                ense = np.unique(xon_match.ensembl_exon_id).astype('|S')[0].decode('utf-8')
                                outls.append('_'.join([str(xon_no), ense, genels[i]]))
                    exonNo1, ense1 = outls[0].split('_')[0:2]
                    exonNo2, ense2 = outls[1].split('_')[0:2]
                    newls = [gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, exonNo1, exonNo2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, warning, confidence, reading_frame, ensg1, ensg2, enst1, enst2, filters, s2_gene, s2_descr, ense1, ense2]
                    print('\t'.join(newls))



