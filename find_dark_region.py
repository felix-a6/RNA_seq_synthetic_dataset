import pysam
from pyfaidx import Fasta
from gtfparse import read_gtf
from multiprocessing import Pool
import argparse
import pickle

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

def find_dark_region(batches, n, mapping_file, gtf_file):
    zero_reads = []
    no_reads = []
    for gene in list(gtf_file.items())[n:batches[batches.index(n)+1]]:
        region =  stat_coverage_refimpl(mapping_file, chrom=gene[1][0], start=gene[1][1], end=gene[1][2], one_based=False)
        region = [x['reads_all'] for x in region]
        if not region:
            no_reads.append(region)
            continue
        if min(region) == 0:
            zero_reads.append(gene[0])
    pickle.dump(zero_reads, open('zero_reads.pkl', 'wb'))

def multi_pross_function(batches):
    for n in batches:
        if batches.index(n) == len(batches)-1:
            break
    find_dark_region(batches, n, mapping_file, gtf_file)    

def get_batches(n_batches, gtf_file):
    batches = [n for n in range(0, len(gtf_file), len(gtf_file)/n_batches)]
    return batches

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate syntetic rna seq dataset')
    parser.add_argument('mapping_file', type=str,
                        help='file containing mapping output')
    parser.add_argument('gtf_file', type=str,
                        help='pickle transcripts infos')

    args = parser.parse_args()
    
    mapping_file = pysam.AlignmentFile(args.mapping_file)
    gtf_file = pickle.load(open(args.gtf_file, 'rb'))

    batches = get_batches(24, gtf_file)
        
    #Multiprocessing
    with Pool(24) as p:
        p.starmap(multi_pross_function, batches)