
# ---------------------------
import sys
from sys import argv
from argparse import ArgumentParser
import os.path
import h5py
import logging

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

# -----------------------------------------------------------------------------
def get_exons_for_transcript(h5_file, transcript_id, position):
    _, pos = position.split(':')
    pos = int(pos)
    with h5py.File(h5_file, 'r') as h5f:
        if transcript_id not in h5f:
            logging.warning(f"Transcript {transcript_id} not found in the HDF5 file.")
            return 'OUT', 'NA'
        transcript_group = h5f[transcript_id]
        exons = []
        for exon_id in transcript_group.keys():
            exon_data = transcript_group[exon_id]
            if isinstance(exon_data, h5py.Dataset):
                exon_values = exon_data[()]
                if len(exon_values) == 3:
                    exon_start, exon_end, exon_number = exon_values
                    exons.append((int(exon_start), int(exon_end), int(exon_number), exon_id))
        if not exons:
            logging.warning(f"No exons found for transcript {transcript_id}")
            return 'NA', 'NA'
        exons.sort(key=lambda x: x[0])  # Sort exons by start position
        for start, end, number, exon_id in exons:
            if start <= pos <= end:
                return str(number), exon_id
        # Check if position is in an intron
        for i in range(len(exons) - 1):
            if exons[i][1] < pos < exons[i+1][0]:
                return 'NA', 'NA'
        # Position is outside the transcript range
        if pos < exons[0][0] or pos > exons[-1][1]:
            return 'OUT', 'NA'
        # This should not happen, but just in case
        return 'NA', 'NA'

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", dest="infus",
                        help="input fusion file name", metavar="<input>")
    parser.add_argument("-d", "--h5db", dest="h5in",
                        help="hdf5 file name", metavar="<h5>")
    args = parser.parse_args()
    fusion_file = args.infus
    h5 = args.h5in

    if not fusion_file:
        print("Please specify a valid fusion file")
        raise SystemExit(1)
    if not h5:
        print("Please specify a valid hdf5 db file name")
        raise SystemExit(1)

    logging.info(f"Processing fusion file: {fusion_file}")
    logging.info(f"Using HDF5 file: {h5}")

    with h5py.File(h5, 'r') as f:
        with open(fusion_file, 'rt') as fus:
            for line_num, line in enumerate(fus, 1):
                line = line.strip('\n')
                if line.startswith('ge'):
                    header = '\t'.join(["gene1", "gene2", "strand1", "strand2", "breakpoint1", "breakpoint2", "exonNo1", "exonNo2", "site1", "site2", "typ", "split_reads1", "split_reads2", "disco_mates", "cov1", "cov2", "warning", "confidence", "reading_frame", "gene_id1", "gene_id2", "transcript_id1", "transcript_id2", "filters", "s2_gene", "s2_descr", "exon_id1", "exon_id2"])
                    print(header)
                else:
                    lspl = line.split('\t')
                    # logging.debug(f"Original split line (length {len(lspl)}): {lspl}")

                    # Ensure we always have 26 fields, filling with empty strings if necessary
                    lspl += [''] * (26 - len(lspl))
                    # logging.debug(f"Padded split line (length {len(lspl)}): {lspl}")

                    [gene1, gene2, strand1, strand2, breakpoint1, breakpoint2, site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, warning, confidence, reading_frame, ensg1, ensg2, enst1, enst2, filters, s2_gene, s2_descr] = lspl[:24]

                    # logging.debug(f"Unpacked values: gene1={gene1}, gene2={gene2}, ..., s2_gene={s2_gene}, s2_descr={s2_descr}")

                    ensts = [enst1, enst2]
                    enstls = [str(x) if x.strip() else 'NA' for x in ensts]
                    breakpoints = [breakpoint1, breakpoint2]

                    # logging.debug(f"ENSTs: {enstls}")
                    # logging.debug(f"Breakpoints: {breakpoints}")

                    outls = []
                    exon_ids = []
                    for i, enst in enumerate(enstls):
                        if enst == 'NA':
                            outls.append('OUT')
                            exon_ids.append('NA')
                        else:
                            if enst in f:
                                exon_number, exon_id = get_exons_for_transcript(args.h5in, enst, breakpoints[i])
                                outls.append(exon_number)
                                exon_ids.append(exon_id)
                            else:
                                logging.warning(f"Transcript {enst} not found in HDF5 file.")
                                outls.append('OUT')
                                exon_ids.append('NA')

                    # logging.debug(f"outls: {outls}")
                    # logging.debug(f"exon_ids: {exon_ids}")

                    out_fields = [gene1, gene2, strand1, strand2, breakpoint1, breakpoint2] + outls + [site1, site2, typ, split_reads1, split_reads2, disco_mates, cov1, cov2, warning, confidence, reading_frame, ensg1, ensg2, enst1, enst2, filters, s2_gene, s2_descr] + exon_ids

                    # logging.debug(f"Final out_fields (length {len(out_fields)}): {out_fields}")

                    if len(out_fields) != 28:
                        logging.warning(f"Incorrect number of output fields on line {line_num}. Expected 28, got {len(out_fields)}")
                        out_fields = out_fields[:28]  # Truncate if too long
                        out_fields += ['NA'] * (28 - len(out_fields))  # Pad if too short

                    out_line = '\t'.join(out_fields)
                    print(out_line)
