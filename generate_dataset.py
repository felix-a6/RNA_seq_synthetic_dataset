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

# module number_reads.py

def nombre_reads_random(bins, numéro_bins):
    debut_bins = bins[numéro_bins]
    fin_bins = bins[numéro_bins+1] 
    return random.randint(debut_bins, fin_bins)  


# module transcript_selection.py

import number_reads

def selection_transcrit(bins, count):
    selected_transcrits_copie = selected_transcrits.copy()
    numéro_bins = 0
    reads_par_transcrit = dict()
    for count in count:
        x=0
        while x<count:
            len_selected_transcrits = len(selected_transcrits_copie)
            transcrit = selected_transcrits_copie.pop(random.randint(0, len_selected_transcrits-1))
            reads_par_transcrit.update({transcrit:nombre_reads_random(bins_made, numéro_bins)})
            x = x+1
        numéro_bins = numéro_bins+1
    return reads_par_transcrit

# module quality_curve_shape.py

def courbe_toutes_positions(scale):
    shape = 2.
    s = np.random.gamma(shape, scale, total_reads)
    count_quality, bins_quality, _ = plt.hist(s, bins = quality_bins_made, density=False)
    return count_quality

# module quality_string.py

def get_quality_string(len_read):
    return ''.join([np.random.choice(quality_levels[:-2], p = pmfs[n]) for n in range(0, len_read)])

# module get_read.py

def make_reads(nom_transcrit, nombre_reads, len_read):
    reads=[]
    transcrit = transcrits_ens[nom_transcrit]
    for _ in range(nombre_reads):
        try:
            position_read = random.randint(0,len(transcrit)-len_read-1)
            reads.append(transcrit[position_read:(position_read+len_read)])
        except Exception as e:
            print(len(transcrit), nom_transcrit)
            raise e
    return reads

# module get_all_reads.py

import get_read
import quality_string

def all_reads(len_read, nom_transcrit_nombre_reads, out_file):
    with open(out_file,'w') as file_fastq:
        for n, (nom_transcrit, nombre_reads) in enumerate(nom_transcrit_nombre_reads):
            reads = make_reads(nom_transcrit, nombre_reads, len_read)
            for read in reads:
                quality_string = get_quality_string(len_read)
                file_fastq.write('\n'.join(['@' + read.fancy_name, str(read), '+', quality_string]) + '\n')

# module helper.py

import get_all_reads

def helper_reads_writer(batch_num):
    all_reads(75, batches[batch_num], f'/mnt/c/Users/xilef/Documents/Stage/Bioinfo/Reads_simulation_{batch_num}.fastq')
