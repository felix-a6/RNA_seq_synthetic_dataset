import random
import numpy as np
import pandas as pd
from math import *
from pyfaidx import Fasta
import itertools as itt
import pickle
import matplotlib.pyplot as plt
import scipy.special as sps
import os
from multiprocessing import Pool
from os.path import join as opj
import argparse
# TODO randomize batches quality strings
ascii_chars = ['!',

 '"',
 '#',
 '$',
 '%',
 '&',
 "'",
 '(',
 ')',
 '*',
 '+',
 ',',
 '-',
 '.',
 '/',
 '0',
 '1',
 '2',
 '3',
 '4',
 '5',
 '6',
 '7',
 '8',
 '9',
 ':',
 ';',
 '<',
 '=',
 '>',
 '?',
 '@',
 'A',
 'B',
 'C',
 'D',
 'E',
 'F',
 'G',
 'H',
 'I',
 'J',
 'K',
 'L',
 'M',
 'N',
 'O',
 'P',
 'Q',
 'R',
 'S',
 'T',
 'U',
 'V',
 'W',
 'X',
 'Y',
 'Z',
 '[',
 '\\',
 ']',
 '^',
 '_',
 '`',
 'a',
 'b',
 'c',
 'd',
 'e',
 'f',
 'g',
 'h',
 'i',
 'j',
 'k',
 'l',
 'm',
 'n',
 'o',
 'p',
 'q',
 'r',
 's',
 't',
 'u',
 'v',
 'w',
 'x',
 'y',
 'z',
 '{',
 '|',
 '}',
 '~']
quality_levels = ascii_chars[:41][::-1]

def nombre_reads_random(bins, numero_bins):
    debut_bins = bins[numero_bins]
    fin_bins = bins[numero_bins+1] 
    return random.randint(debut_bins, fin_bins)  

def generate_reads_distribution(output_dir, max_n_reads = 4000, n_bins = 251, n_transcript = 40000, high = (1000, 500), low = (0, 500)):
    """combine 2 normal distributions high = (mean, standard dev) low = (mean, standard dev)
    """
    #Bins used for model
    bins_made = np.linspace(0, max_n_reads, n_bins)
    
    #Poids des transcrits Distribution normale
    s_high = np.random.normal(high[0], high[1], n_transcript) 
    s_low = np.random.normal(low[0], low[1], n_transcript) 
    count_high, bins_high, _ = plt.hist(s_high, bins = bins_made, density=False)
    count_low, bins_low, _ = plt.hist(s_low, bins = bins_made, density=False)

    #Sélection du plus grand compte des deux courbes
    count_s = []
    for n in range(n_bins - 1):
        count_s.append(max([count_low[n], count_high[n]]))
    count_s.append(0)
    plt.plot(bins_made, count_s)
    plt.savefig(opj(output_dir, 'model_read_distribution.png'))
    plt.show()

    #Compte ajusté par pmf
    count_sum = sum(count_s)
    count_adjusted = [(x/count_sum)*n_transcript for x in count_s]
    
    return bins_made, count_adjusted

def selection_transcrit(bins_made, counts, selected_transcrits): # TODO Better function name
    selected_transcrits_copie = selected_transcrits.copy()
    reads_par_transcrit = dict()
    for numero_bins, count in enumerate(counts):
        x = 0
        while x<count:
            len_selected_transcrits = len(selected_transcrits_copie)
            transcrit = selected_transcrits_copie.pop(random.randint(0, len_selected_transcrits-1))
            reads_par_transcrit.update({transcrit:nombre_reads_random(bins_made, numero_bins)})
            x = x+1
    return reads_par_transcrit

def get_quality_pmfs(output_dir, total_reads, len_read = 75, quality_bins_made = np.linspace(0, 40, 40)):
    positions = list(range(1, len_read+1))
    scales = [0.5 + position * (2/len_read) for position in positions]
    data = []
    pmfs = []
    
    def courbe_toutes_positions(scale):
        shape = 2.
        s = np.random.gamma(shape, scale, total_reads)
        count_quality, bins_quality, _ = plt.hist(s, bins = quality_bins_made, density=False)
        return count_quality
    
    for scale in scales:
        count = courbe_toutes_positions(scale)
        data.append(count)
        pmfs.append(count/sum(count))
    pickle.dump(data, open(opj(output_dir, 'read_quality_distribution.pkl'), 'wb'))
    return pmfs

def get_quality_string(len_read, pmfs):
    return ''.join([np.random.choice(quality_levels[:-2], p = pmfs[n]) for n in range(0, len_read)])

def make_reads(nom_transcrit, nombre_reads, len_read, transcrits_sequences):
    reads=[]
    transcrit = str(transcrits_sequences[nom_transcrit])
    for _ in range(nombre_reads):
        position_read = random.randint(0,len(transcrit)-len_read-1)
        read = transcrit[position_read:(position_read+len_read)]
        if len(read) == 75:
            reads.append((read, position_read))
        else:
            print(f'len_err : {nom_transcrit}, len_read:{len(read)}, position_read:{position_read}, nombre_reads:{nombre_reads}', flush = True)
    return reads

def get_batches(n_batches, transcript_reads, output_dir, transcrits_sequences, pmfs, len_reads):
    batch_size = int(np.ceil(1/n_batches * len(transcript_reads)))
    transcript_reads_shuffle = list(transcript_reads.items())
    random.shuffle(transcript_reads_shuffle) #TODO enlever list comprehention
    batches = []
    for n in range(0,  len(transcript_reads_shuffle), batch_size):
        out_file = opj(output_dir, f'Reads_simulation_{n}.fastq')
        batch = transcript_reads_shuffle[n:n+batch_size]
        batches.append((len_reads, batch, out_file, transcrits_sequences, pmfs))  
    return batches

def write_batch_reads(len_read, batch, out_file, transcrits_sequences, pmfs):
    with open(out_file,'w') as file_fastq:
        print_counter = 0
        for n, (nom_transcrit, nombre_reads) in enumerate(batch):
            if print_counter == 100:
                print(nom_transcrit, nombre_reads, flush = True)
            reads = make_reads(nom_transcrit, nombre_reads, len_read, transcrits_sequences)
            for read, position_read in reads:
                quality_string = get_quality_string(len_read, pmfs)
                read_name = f'@{nom_transcrit}:{position_read}'
                file_fastq.write('\n'.join([read_name, read, '+', quality_string]) + '\n')
                if print_counter == 100:
                    print(read_name, flush = True)
                    print_counter = 0
            print_counter += 1

def combine_files(output_dir):
    with open(opj(output_dir, 'Reads_simulation.fastq'), 'w') as f_write:
        for file in os.listdir(output_dir):
            if 'Reads_simulation_' in file:
                with open(opj(output_dir, file), 'r') as f_read:
                    for line in f_read:
                        f_write.write(line)
                # TODO delete batch file

def truncate(seq, length=80):
    return '\n'.join([seq[n:n+length] for n in range(0, len(seq), length)])

def get_dictionnary(transcrits_sequences):
    dict_transcrits_sequences = dict()
    for transcrit, sequence in transcrits_sequences.items():
        dict_transcrits_sequences.update({transcrit:str(sequence)})
    return dict_transcrits_sequences

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate syntetic rna seq dataset')
    parser.add_argument('fasta', type=str,
                        help='fasta containing all transcript sequences')
    parser.add_argument('output_dir', type=str,
                        help='directory for all outputs')
    parser.add_argument('n_batches', type=int,
                        help='number of batches and CPUs')


    args = parser.parse_args()
     
    transcrits_sequences = Fasta(args.fasta)
    output_dir = args.output_dir
    n_batches = args.n_batches
    selected_transcrits = list(transcrits_sequences.keys())
    bins_made, count_adjusted = generate_reads_distribution(output_dir)

    if not os.path.exists(opj(output_dir, 'transcript_reads.pkl')):
        transcript_reads=selection_transcrit(bins_made, count_adjusted, selected_transcrits)
        pickle.dump(transcript_reads, open(opj(output_dir, 'transcript_reads.pkl'), 'wb'))
    else:
        transcript_reads = pickle.load(open(opj(output_dir, 'transcript_reads.pkl'), 'rb'))

    total_reads = sum(transcript_reads.values())

    if os.path.exists(opj(output_dir, 'quality_pmfs.pkl')):
        pmfs = pickle.load(open(opj(output_dir, 'quality_pmfs.pkl'), 'rb'))
    else:
        pmfs =get_quality_pmfs(output_dir, total_reads)
        pickle.dump(pmfs, open(opj(output_dir, 'quality_pmfs.pkl'), 'wb'))

    transcrits_sequences = get_dictionnary(transcrits_sequences)

    batches = get_batches(n_batches, transcript_reads, output_dir, transcrits_sequences, pmfs, 75)

    #Multiprocessing
    with Pool(n_batches) as p:
        p.starmap(write_batch_reads, batches)
    combine_files(output_dir)