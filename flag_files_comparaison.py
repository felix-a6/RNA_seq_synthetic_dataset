import pickle
from matplotlib_venn import venn3
import argparse

def create_set(pickle_file):
    multimapped_reads =  [read.split('\t')[0] for read in pickle_file]
    return set(multimapped_reads)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate syntetic rna seq dataset')
    parser.add_argument('pickle_150', type=str,
                        help='pickle file for reads 150')
    parser.add_argument('pickle_100', type=str,
                        help='pickle file for reads 100')
    parser.add_argument('pickle_75', type=int,
                        help='pickle file for reads 75')
    parser.add_argument('output_dir', type=int,
                        help='output directory')

    output_dir = args.output_dir

    set_150 = create_set(args.pickle_150)
    set_100 = create_set(args.pickle_100)
    set_75 = create_set(args.pickle_75)

    venn_diagram = venn3([set_75, set_100, set_150], set_labels = ['75', '100', '150'])

    pickle.dump(venn_diagram, open(f'{output_dir}venn_diagram.pkl', 'wb'))