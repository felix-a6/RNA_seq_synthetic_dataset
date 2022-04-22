from pyfaidx import Fasta
import pickle
from gtfparse import read_gtf
import pysam
import os
import numpy as np
import argparse

def get_1_read_count(input_path, files_name):
    all_1_read_count = []
    for file in os.listdir(input_path):
        if files_name not in file: continue
        zero_reads = pickle.load(open(f'{input_path}{file}', 'rb'))
        for trans in zero_reads:
            all_1_read_count.append(trans)
    return all_1_read_count

def get_dark_genes_coords(gtf_file, input_path, files_name):
    chroms = set([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT'])
    gene_coordinates = gtf_file[(gtf_file['feature']=='gene') & (gtf_file['seqname'].isin(chroms))][['seqname','strand', 'start', 'end', 'gene_id']].to_dict(orient='records')
    dark_gene_coords = []
    for coords in gene_coordinates:
        if coords['gene_id'] in get_1_read_count(input_path, files_name):
            dark_gene_coords.append(coords)
    return dark_gene_coords

def normalise_coords(one_based, start, end):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    return start, end

def stat_coverage_refimpl(samfile, chrom=None, start=None, end=None,
                          one_based=True):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = col.pileups
        yield {'chrom': chrom, 'pos': pos, 'reads_all': len(reads)}

def get_dark_edges(gene_coords, pos):
    dark_edges = []
    if pos[0] > gene_coords['start']:
        dark_edges.append((gene_coords['start'], pos[0]))
    if pos[len(pos)-1] < gene_coords['end']:
        dark_edges.append((pos[len(pos)-1],  gene_coords['end']))
    return dark_edges

def get_dark_region_coords(gene_coords, mapping_file):
    region =  stat_coverage_refimpl(mapping_file, chrom=gene_coords['seqname'], start=gene_coords['start'], end=gene_coords['end'], one_based=False)
    region = list(region)

    if not region:
        print('no_reads:', gene_coords)
        return []

    (pos, read_counts) = zip(*[(nt['pos'], nt['reads_all']) for nt in region])

    pos = np.array(pos)
    gaps = pos[1:]-pos[:-1]>1
    start = pos[:-1][np.where(gaps)]
    end = pos[1:][np.where(gaps)]
    dark_coords = list(zip(start, end))
    
    for edge in get_dark_edges(gene_coords, pos):
        dark_coords.append(edge)
    
    return dark_coords

def get_all_dark_region(gtf_file, input_path, files_name, mapping_file, output_dir):
    all_dark_regions = dict()
    for gene_coords in get_dark_genes_coords(gtf_file, input_path, files_name):
        all_dark_regions.update({gene_coords['gene_id']:(gene_coords['seqname'], gene_coords['strand'], get_dark_region_coords(gene_coords, mapping_file))})
    pickle.dump(all_dark_regions, open(f'{output_dir}dark_region_coords_new.pkl', 'wb'))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Dark regions data analysis')
    parser.add_argument('gtf_file', type=str,
                        help='path to gtf_file')
    parser.add_argument('output_dir', type=str,
                        help='output direction')
    parser.add_argument('input_path', type=str,
                        help='path to input zeros reads')
    parser.add_argument('files_name', type=str,
                        help='commun name zeros reads')
    parser.add_argument('mapping_file', type=str,
                        help='path to mapping file')

    args = parser.parse_args()
    
    gtf_file = read_gtf(args.gtf_file)
    output_dir = args.output_dir
    input_path = args.input_path
    files_name = args.files_name
    mapping_file = pysam.AlignmentFile(args.mapping_file, 'rb')


    get_all_dark_region(gtf_file, input_path, files_name, mapping_file, output_dir)