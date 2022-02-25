import pysam
from matplotlib import pyplot as plt
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
    chroms = set([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT'])
    gene_coordinates = gtf_file[(gtf_file['feature']=='gene') & (gtf_file['seqname'].isin(chroms))][['seqname', 'start', 'end', 'gene_id']].to_dict(orient='records')
    zero_reads = []
    no_reads = []
    for gene in gene_coordinates[n:batches[batches.index(n)+1]]:
        region =  stat_coverage_refimpl(mapping_file, chrom=gene['seqname'], start=gene['start'], end=gene['end'], one_based=False)
        region = [x['reads_all'] for x in region]
        if not region:
            no_reads.append(region)
            continue
        if min(region) == 0:
            zero_reads.append(gene['gene_id'])
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
                        help='file containing transcripts infos')

    args = parser.parse_args()
    
    mapping_file = pysam.AlignmentFile(args.mapping_file)
    gtf_file = read_gtf(args.gtf_file)

    batches = get_batches(24, gtf_file)
        
    #Multiprocessing
    with Pool(24) as p:
        p.starmap(multi_pross_function, batches)