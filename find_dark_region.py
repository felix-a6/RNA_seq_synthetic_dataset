import pysam
from pyfaidx import Fasta
from multiprocessing import Pool
import argparse
import pickle
import numpy as np

def normalise_coords(one_based, start, end):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    return start, end

def stat_coverage_refimpl(samfile, chrom=None, start=None, end=None,
                          one_based=False):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = col.pileups
        yield {'chrom': chrom, 'pos': pos, 'reads_all': len(reads)}

def find_dark_region(batch, batch_num, mapping_file, output_dir):
    zero_reads = []
    no_reads = []
    for gene_name, coords in batch:
        chrom, start, end = coords
        region =  stat_coverage_refimpl(mapping_file, chrom, start, end, one_based=False)
        region = [x['reads_all'] for x in region]
        if not region:
            no_reads.append(region)
            continue
        if min(region) == 0:
            zero_reads.append(gene_name)
    pickle.dump(zero_reads, open(f'{output_dir}/zero_reads{batch_num}.pkl', 'wb'))

def get_batches(n_batches, gtf_file, mapping_file, output_dir):
    batch_size = int(np.ceil(len(gtf_file)/n_batches))
    batches = [
        (list(gtf_file.items())[n:n+batch_size], n, mapping_file, output_dir)
        for n in range(0, len(gtf_file), batch_size)
        ]
    return batches

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate syntetic rna seq dataset')
    parser.add_argument('mapping_file', type=str,
                        help='file containing mapping output')
    parser.add_argument('gtf_file', type=str,
                        help='pickle transcripts infos')
    parser.add_argument('output_dir', type=str,
                        help='output direction')

    args = parser.parse_args()
    
    mapping_file = pysam.AlignmentFile(args.mapping_file, 'rb')
    gtf_file = pickle.load(open(args.gtf_file, 'rb'))

    batches = get_batches(24, gtf_file, mapping_file, args.output_dir)
        
    #Multiprocessing
    with Pool(24) as p:
        p.starmap(find_dark_region, batches)