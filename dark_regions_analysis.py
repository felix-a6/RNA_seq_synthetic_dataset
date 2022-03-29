import pickle
import os
import pysam
from gtfparse import read_gtf
from matplotlib import pyplot as plt
from pyfaidx import Fasta
from multiprocessing import Pool
import numpy as np
from sympy import Interval, Union
import matplotlib as mpl
from os.path import join as opj
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
    gene_coordinates = gtf_file[(gtf_file['feature']=='gene') & (gtf_file['seqname'].isin(chroms))][['seqname', 'start', 'end', 'gene_id']].to_dict(orient='records')
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

    pos, read_counts = zip(*[(nt['pos'], nt['reads_all']) for nt in region])

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
        all_dark_regions.update({gene_coords['gene_id']:(gene_coords['seqname'], get_dark_region_coords(gene_coords, mapping_file))})
    pickle.dump(all_dark_regions, open(f'{output_dir}dark_region_coords.pkl', 'wb'))

def get_exons_coords_dict(gtf_file, output_dir):
    chroms = set([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT'])
    exons_coordinates = gtf_file[(gtf_file['feature']=='exon') & (gtf_file['seqname'].isin(chroms))][['seqname', 'start', 'end', 'gene_id']].to_dict(orient='records')
    gene_exons_dict = dict()
    for exon in exons_coordinates:
        if exon['gene_id'] in (pickle.load(open(f'{output_dir}dark_region_coords', 'rb'))).keys():
            if exon['gene_id'] not in gene_exons_dict:
                gene_exons_dict.update({exon['gene_id']:[(exon['start'],exon['end'])]})
            if exon['gene_id'] in gene_exons_dict:
                gene_exons_dict[exon['gene_id']].append((exon['start'],exon['end']))
    return gene_exons_dict
            
def get_introns_dict(gtf_file, output_dir):
    introns_dict = dict()
    gene_exons_dict = get_exons_coords_dict(gtf_file, output_dir)
    for exons in gene_exons_dict:
        gene_start, gene_end = min([x[0] for x in gene_exons_dict[exons]]), max([x[1] for x in gene_exons_dict[exons]])
        gene_length = gene_end-gene_start
        gene_map = np.zeros(gene_length)
        for start, end in gene_exons_dict[exons]:
            gene_map[start-gene_start:end-gene_start] = 1

        introns_map = np.where(gene_map)[0]
        gaps = introns_map[1:]-introns_map[:-1]>1
        start = introns_map[:-1][np.where(gaps)]
        end = introns_map[1:][np.where(gaps)]
        introns_coords = list(zip(start+1+gene_start, end-1+gene_start))
        introns_dict.update({exons:introns_coords})
    return introns_dict

def remove_introns(dark_start, dark_end, introns_list):
    introns_in_dark_length = 0
    for start_intron, end_intron in introns_list:
        if start_intron in range(dark_start, dark_end) and end_intron in range(dark_start, dark_end):
            introns_in_dark_length = (end_intron - start_intron) + introns_in_dark_length
    return introns_in_dark_length

def get_dark_region_length(gtf_file, output_dir):
    introns_dict = get_introns_dict(gtf_file, output_dir)
    dark_region_coords = pickle.load(open(f'{output_dir}dark_region_coords', 'rb'))
    dark_region_length = []
    for gene, coords in dark_region_coords.items():
        chrom, dark_regions = coords
        for dark_region in dark_regions:
            start, end = dark_region
            introns_in_dark_length = remove_introns(start, end, introns_dict[gene])
            dark_region_length.append((gene, end - start - introns_in_dark_length))
    return dark_region_length

def get_length_hist(gtf_file, output_dir):
    len_dark = []
    for dark_region in get_dark_region_length(gtf_file, output_dir):
        _, lenght = dark_region
        len_dark.append(lenght)
    plt.hist(np.log10(np.array(len_dark)),  bins = 50)
    plt.title('dark_regions_length')
    plt.xlabel('Log10_lenght')
    plt.ylabel('count')
    plt.gca().set_yscale('log')
    plt.savefig(opj(output_dir, 'dark_region_length_hist.png'))
    plt.show()

def get_biotype_dict(biotype_file, output_dir):
    with open(biotype_file, 'r') as biotype:
        dark_region_coords = pickle.load(open(f'{output_dir}dark_region_coords', 'rb'))
        biotype_dict = dict()
        biotype = biotype.readlines()
        for gene, coords in dark_region_coords.items():
            chrom, dark_region = coords
            if dark_region == []:
                continue
            for line in biotype:
                if gene in line:
                    line = '\n'.join(line.split('\t',)).split('\n',)
                    if gene not in biotype_dict:
                        biotype_dict.update({gene:[line[3]]})
                    if gene in biotype_dict:
                        if line[3] not in biotype_dict[gene]:
                            biotype_dict[gene].append(line[3])
        return biotype_dict

def get_all_gene_biotype(biotype_file, output_dir):
    biotype_groups = {
    'processed_transcript':'non_coding',
    'lincRNA':'non_coding',
    'antisense':'non_coding',
    'retained_intron':'non_coding',
    'sense_intronic':'non_coding',
    'sense_overlapping':'non_coding',
    
    'transcribed_unprocessed_pseudogene':'pseudogene',
    'unprocessed_pseudogene':'pseudogene',
    'processed_pseudogene':'pseudogene',
    'transcribed_processed_pseudogene':'pseudogene',
    'pseudogene':'pseudogene',
    'polymorphic_pseudogene':'pseudogene',
    'transcribed_unitary_pseudogene':'pseudogene',
    'rRNA_pseudogene':'pseudogene',
    
    'protein_coding':'protein_coding',
    
    'nonsense_mediated_decay':'non_coding',
    'miRNA':'non_coding',
    'TEC':'non_coding',
    'rRNA':'non_coding',
    'bidirectional_promoter_lncRNA':'non_coding',
    'snRNA':'non_coding',
    'non_stop_decay':'non_coding',
    'snoRNA':'non_coding'
}
    gene_biotype_dict = dict()
    for gene, biotypes in (get_biotype_dict(biotype_file, output_dir)).items():
        biotypes_groups =[]
        for biotype in biotypes:
            if biotype_groups[biotype] not in biotypes_groups:
                biotypes_groups.append(biotype_groups[biotype])
        gene_biotype_dict.update({gene:biotypes_groups})
    return gene_biotype_dict

def dark_region_by_biotypes(biotype_file, output_dir):
    protein_coding_group = []
    pseudogene_group = []
    protein_coding_and_pseudogene_group = []
    non_coding_group = []
    gene_biotype_dict = get_all_gene_biotype(biotype_file, output_dir)
    for gene in gene_biotype_dict:
        if 'protein_coding' in gene_biotype_dict[gene] and 'pseudogene' in gene_biotype_dict[gene]:
            protein_coding_and_pseudogene_group.append(gene)
        elif 'pseudogene' in gene_biotype_dict[gene]:
            pseudogene_group.append(gene)
        elif 'protein_coding' in gene_biotype_dict[gene]:
            protein_coding_group.append(gene)
        else:
            non_coding_group.append(gene)
    return protein_coding_group, pseudogene_group, protein_coding_and_pseudogene_group, non_coding_group

def get_length_by_biotypes(biotype_list, output_dir):
    dark_region_length = []
    dark_region_coords = pickle.load(open(f'{output_dir}dark_region_coords', 'rb'))
    for gene in biotype_list:
        chrom, dark_regions = dark_region_coords[gene]
        for dark_region in dark_regions:
            start, end = dark_region
            dark_region_length.append(end - start)
    return dark_region_length

def get_figs_biotypes(biotype_file, output_dir):
    protein_coding_group, pseudogene_group, protein_coding_and_pseudogene_group, non_coding_group = dark_region_by_biotypes(biotype_file, output_dir)
    distribution = [len(protein_coding_group), len(pseudogene_group), len(protein_coding_and_pseudogene_group), len(non_coding_group)]
    interesting_groups = ['protein_coding', 'pseudogene', 'both_biotypes', 'non_coding_group']
    y_pos = np.arange(len(interesting_groups))

    plt.bar(y_pos, distribution, align='center', alpha=0.5)
    plt.xticks(y_pos, interesting_groups)
    plt.ylabel('count')
    plt.title('genes biotypes containing dark regions')
    plt.savefig(opj(output_dir, 'biotypes_distribution.png'))
    plt.show()

    protein_coding_length = get_length_by_biotypes(protein_coding_group, output_dir)
    pseudogene_length = get_length_by_biotypes(pseudogene_group, output_dir)
    non_coding_length = get_length_by_biotypes(non_coding_group, output_dir)

    plt.hist(np.log10(np.array(protein_coding_length)), bins = 50, histtype = 'step', label = 'protein_coding')
    plt.hist(np.log10(np.array(pseudogene_length)), bins = 50, histtype = 'step', label = 'pseudogene')
    plt.hist(np.log10(np.array(non_coding_length)), bins = 50, histtype = 'step', label = 'non_coding')
    plt.gca().set_yscale('log')
    plt.legend()
    plt.savefig(opj(output_dir, 'length_by_biotypes.png'))
    plt.show()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Dark regions data analysis')
    parser.add_argument('gtf_file', type=str,
                        help='gtf file path')
    parser.add_argument('input_path', type=str,
                        help='path to dark regions files')
    parser.add_argument('files_name', type=str,
                        help='commun name of dark regions files')
    parser.add_argument('mapping_file', type=str,
                        help='output from mapping, aligned by coords')
    parser.add_argument('output_dir', type=str,
                        help='output direction')
    parser.add_argument('biotype_file', type=str,
                        help='txt for transcripts biotypes')

    args = parser.parse_args()
    
    gtf_file = read_gtf(args.gtf_file)
    input_path = args.input_path
    files_name = args.files_name
    mapping_file = args.mapping_file
    output_dir = args.output_dir
    biotype_file = args.biotype_file

    get_all_dark_region(gtf_file, input_path, files_name, mapping_file, output_dir)
    get_length_hist(gtf_file, output_dir)
    get_figs_biotypes(biotype_file, output_dir)