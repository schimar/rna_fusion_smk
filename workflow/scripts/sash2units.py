#! /usr/bin/python


from sys import argv
from itertools import islice, dropwhile
from glob import glob
import re
# -----------------------------

sampleSheet = argv[1]
fqpath = argv[2]

# --------

def get_sid_fullpath(fqpath, sid):
    for fq in glob(f"{fqpath}*fastq.gz"):
        p = re.compile(f"{fqpath}{sid}[\w]+R1[\w.]+")
        if re.search(p, fq) is None:
            continue
        else:
            matchfq = re.search(p, fq).group(0)
            if re.search("L001", matchfq) is not None:
                unit = "L001"
            elif re.search("L002", matchfq) is not None:
                unit = "L002"
            else:
                continue
            return matchfq, unit

# -----------------------------------------------------------------

if __name__ == "__main__":

    with open(sampleSheet, 'rt') as sas:
        itskip = islice(sas, 20, None)
        for line in itskip:
            lspl = line.strip('\n').split(',')
            if lspl[0] == "Sample_ID":
                header = '\t'.join(lspl)
                #print(header)
            else:
                sid = lspl[0]
                adapters = ' '.join([lspl[6], lspl[8]])
                fqR1, unit = get_sid_fullpath(fqpath, sid)
                fqR2 = fqR1.replace('_R1_', '_R2_')
                print(f"{sid}\t{unit}\t{adapters}\t{fqR1}\t{fqR2}")




#for _ in range(17):
#        next(f)
#    for line in f:
#        # do stuff


