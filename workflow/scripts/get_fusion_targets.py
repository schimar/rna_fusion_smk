
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

# Usage: scripts/get_fusion_targets.py -i fusions.tsv -t targets.tsv

# -----------------------------------------------------------------------------


# ---------------------------------- #

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="infus",
                        help="input fusion file name", metavar="<input>")
    parser.add_argument("-t", "--targets", dest="intargets",
                        help="input target file name", metavar="<input>")

    #parser.add_argument("-q", "--quiet",
    #                    action="store_false", dest="verbose", default=True,
    #                    help="don't print status messages to stdout")

    args = parser.parse_args()

    fusion_file = args.infus
    target_file = args.intargets

    if not fusion_file: #.exists():
        print("Please specify a valid fusion file")
        raise SystemExit(1)


    targetfusdict = dict()
    with open(target_file, 'rt') as targs:
        for line in targs:
            line = line.strip('\n')
            lspl = line.split('\t')
            fusion = lspl[0]    #.split(':')
            if len(lspl) > 1:
                xtra = lspl[1].split(',')
            else:
                xtra = ['']
            if fusion not in targetfusdict:
                targetfusdict[fusion] = xtra
            else:
                continue



    arrfusdict = dict()
    with open(fusion_file, 'rt') as fus:
        lines = fus.readlines()
        for line in lines:
            line = line.strip('\n')
            if line[0:2] == '#g':
                continue
                #print(line)
            else:
                lspl = line.split()
                #reading_frame = lspl[15]
                gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, confidence, reading_frame, tags, retained_protein_domains, closest_genomic_breakpoint1, closest_genomic_breakpoint2, gene_id1, gene_id2, transcript_id1, transcript_id2, direction1, direction2, filters, fusion_transcript, peptide_sequence, read_id = lspl[0:30]
                fusion = ':'.join([gene1, gene2])
                if not fusion in arrfusdict:
                    arrfusdict[fusion] = [[breakpoint1, breakpoint2]]
                else:
                    arrfusdict[fusion].append([breakpoint1, breakpoint2])

    anotb = list()
    bnota = list()
    ainb = list()
    for key, value in targetfusdict.items():
        if len(value) > 1:
            if len(value) == len(arrfusdict[key]):
                ainb.append(f"{key} (2trans.)")
                #print(f"{key} - {len(value)} fusion transcripts found")
        else:
            if not key in arrfusdict:
                anotb.append(key)
                #print(f"{key} not found in arriba output")
            elif key in arrfusdict:
                ainb.append(key)
                #print(f"{key} fusion present")
            else:
                print("somethings's up")
            #print(key, len(value))

    for key, value in arrfusdict.items():
        if key not in targetfusdict:
            bnota.append(key)
            #print(f"{key} not in list of target fusions")

        #check if list is nested
        # any(isinstance(i, list) for i in unlist(value)))

    if len(anotb) == 0:
        anotb= "0"
    if len(bnota) == 0:
        bnota= "0"

    nl = '\n'
    print(f"targets found in arriba output: {nl}{nl.join(ainb)}", '\n')
    print(f"fusions not in target file: {nl}{nl.join(bnota)}", '\n')
    print(f"targets not in arriba output: {nl}{nl.join(anotb)}", '\n')


    fus.close()
    targs.close()


