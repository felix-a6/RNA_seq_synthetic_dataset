import pickle
import argparse
from collections import defaultdict

def get_chrom_seq(chrom_doc_path):
    with open(chrom_doc_path, 'r') as chrom_doc:
        lines = [line for line in chrom_doc if '>' not in line]
        chrom_seq = ''.join(''.join(lines).split('\n'))
        return chrom_seq

def get_chroms_seq_dict(chroms_files_path, output_dir):
    chroms_seq_dict = dict()
    chroms = set([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT'])
    for chrom in chroms:
        chrom_seq = get_chrom_seq(f'{chroms_files_path}Homo_sapiens.GRCh38.dna.chromosome.{chrom}.fa')
        chroms_seq_dict.update({chrom:chrom_seq})
    pickle.dump(chroms_seq_dict, open(f'{output_dir}chroms_seq_dict.pkl', 'wb'))

def get_dark_region_seq(starts_list, ends_list, chrom_seq):
    darks_exons = []
    for n, (start) in enumerate(starts_list):
        darks_exons.append(chrom_seq[start-140:ends_list[n]+140])
    return ''.join(darks_exons)

def get_coords_list(dark_region_coords, gene_introns): #dark_region_coords 1 dark_region and introns_dict for 1 gene
    (dark_start, dark_end) = dark_region_coords
    starts_list = [dark_start]
    ends_list = []
    for (intron_start, intron_end) in gene_introns:
        if intron_start in range(dark_start-140, dark_end+140) and intron_end in range(dark_start-140, dark_end+140):
            starts_list.append(intron_end)
            ends_list.append(intron_start)
    ends_list.append(dark_end)
    return starts_list, ends_list

def update_seq_dict(gene_name, dark_region_coords, gene_introns, chrom_seq, dark_regions_seq_dict):
    (starts_list, ends_list) = get_coords_list(dark_region_coords, gene_introns)
    dark_region_seq = get_dark_region_seq(starts_list, ends_list, chrom_seq)
    dark_regions_seq_dict.update({f'{gene_name}-{dark_region_coords}':dark_region_seq})
    return dark_regions_seq_dict

def get_all_seq_dict(dark_region_coords_dict, intron_dict, output_dir):
    chroms_seq_dict = pickle.load(open(f'{output_dir}chroms_seq_dict.pkl', 'rb'))
    dark_regions_seq_dict = dict()
    for gene_name in dark_region_coords_dict.keys():
        (chrom_num, dark_coords_list) = dark_region_coords_dict[gene_name]
        gene_introns = intron_dict[gene_name]
        chrom_seq = chroms_seq_dict[chrom_num]
        for dark_coords in dark_coords_list:
            dark_regions_seq_dict = update_seq_dict(gene_name, dark_coords, gene_introns, chrom_seq, dark_regions_seq_dict)

    sorted_dark_region_seq_dict = defaultdict(list)
    for key, val in sorted(dark_regions_seq_dict.items()):
        sorted_dark_region_seq_dict[val].append(key)

    pickle.dump(sorted_dark_region_seq_dict, open(f'{output_dir}sorted_dark_seq_dict.pkl', 'wb'))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Dark regions data analysis')
    parser.add_argument('chroms_files_path', type=str,
                        help='path to chroms files')
    parser.add_argument('output_dir', type=str,
                        help='output direction')
    parser.add_argument('dark_region_coords_dict', type=str,
                        help='path to dark regions coords pickle')
    parser.add_argument('intron_dict', type=str,
                        help='path to intron_dict pickle')

    args = parser.parse_args()
    
    chroms_files_path = args.chroms_files_path
    output_dir = args.output_dir
    dark_region_coords_dict = pickle.load(open(args.dark_region_coords_dict, 'rb'))
    intron_dict = pickle.load(open(args.intron_dict, 'rb'))

    get_chroms_seq_dict(chroms_files_path, output_dir)
    get_all_seq_dict(dark_region_coords_dict, intron_dict, output_dir)