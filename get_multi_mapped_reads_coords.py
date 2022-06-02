import argparse
import pysam
import pickle

def get_coords(bam_file, output_dir):
    coords_list = [(read.qname, read.reference_name, read.reference_start, read.reference_end, read.rlen, read.cigarstring) for read in bam_file.fetch()]
    pickle.dump(coords_list, open(f'{output_dir}multi_mapped_coords.pkl', 'wb'))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Multi mapped coords')
    parser.add_argument('bam_file', type=str,
                        help='path to bam_file')
    parser.add_argument('output_dir', type=str,
                        help='output_directory')

    args = parser.parse_args()

    bam_file = pysam.AlignmentFile(args.bam_file, 'rb')
    output_dir = args.output_dir

    get_coords(bam_file, output_dir)