

from argparse import ArgumentParser
import h5py

def write2hdf5(h5file, attr_dict, start, end):
    transcript_id = attr_dict['transcript_id']
    exon_id = attr_dict['exon_id']
    exon_number = attr_dict.get('exon_number', '0')

    if transcript_id not in h5file:
        h5file.create_group(transcript_id)
    if exon_id not in h5file[transcript_id]:
        exon_data = h5file[transcript_id].create_dataset(exon_id, (3,), dtype='i4')
    else:
        exon_data = h5file[transcript_id][exon_id]

    exon_data[0] = int(start)
    exon_data[1] = int(end)
    exon_data[2] = exon_number

def gtf2exon_h5(gtf, h5):
    with h5py.File(h5, 'w') as h5file:
        with open(gtf, 'rt') as anno:
            transcript0 = None
            for line in anno:
                if line[0] != '#':
                    lspl = line.strip('\n').split('\t')
                    seqname, source, feature, start, end, score, strand, frame, attr = lspl
                    attrspl = attr.split(';')
                    attr_dict = {key_value.strip().split(' ')[0]: key_value.strip().split(' ')[1].strip('"')
                                 for key_value in attrspl if key_value.strip()}
                    if feature == 'exon' and 'transcript_id' in attr_dict and 'exon_id' in attr_dict:
                        if transcript0 != attr_dict['transcript_id']:
                            write2hdf5(h5file, attr_dict, start, end)
                            transcript0 = attr_dict['transcript_id']
                        else:
                            write2hdf5(h5file, attr_dict, start, end)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-g", "--gtf", dest="gtf",
                        help="input gtf file")
    parser.add_argument("-h5", "--hdf5", dest="h5",
                        help="output hdf5 file")
    args = parser.parse_args()

    if not args.gtf:
        print("Please specify a valid GTF file")
        raise SystemExit(1)
    if not args.h5:
        print("Please specify a valid HDF5 output file name")
        raise SystemExit(1)

    gtf2exon_h5(args.gtf, args.h5)
