from pyfaidx import Fasta
import argparse
import pickle

def truncate(seq, length=80):
    return '\n'.join([seq[n:n+length] for n in range(0, len(seq), length)])

def write_all_reads_fasta(transcript_fasta, len_read, out_file, step_size):
    with open(out_file,'w') as all_reads_fasta:
        too_short = []
        for name, seq in transcript_fasta.items():
            if len(seq) < 150 :
                too_short.append(name)
                continue
            reads = [seq[position_read:(position_read+len_read)] for position_read in range(0, len(seq)-len_read, step_size)]
            for read in reads:
                all_reads_fasta.write('>' + name + '\n' + truncate(str(read),length=80) + '\n')
        pickle.dump(too_short, open(f'{out_file}_short', 'wb'))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate syntetic rna seq dataset full coverage')
    parser.add_argument('fasta', type=str,
                        help='fasta containing all transcript sequences')
    parser.add_argument('out_file', type=str,
                        help='fasta where output is')
    parser.add_argument('len_read', type=int,
                        help='length reads')
    parser.add_argument('step_size', type=int, default = 1,
                        help='nummber of nucleotide between reads')

    args = parser.parse_args()

    transcript_fasta = Fasta(args.fasta)

    write_all_reads_fasta(transcript_fasta, args.len_read, args.out_file, args.step_size)