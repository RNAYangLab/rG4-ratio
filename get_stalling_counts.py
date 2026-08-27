import time
import argparse
import pysam
from Bio import SeqIO
from multiprocessing import Pool, cpu_count


def stalling_counting(transcript):
    seq = str(seq_dict[transcript].seq)
    length = len(seq)
    result_dict = {'seq': seq, 'counts': [0] * length, 'coverage': [0] * length}
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(transcript):
            if (read.is_read1 and read.mate_is_mapped and read.is_forward and
                    read.mate_is_reverse and read.reference_name == read.next_reference_name):
                mapping_start = read.reference_start
                result_dict['counts'][mapping_start] += 1
                for pos in read.positions:
                    result_dict['coverage'][pos] += 1
    return (transcript, result_dict)


def extract_transcript_and_counting():
    transcript_list = list(seq_dict.keys())
    with Pool(args.threads) as pool:
        transcript_dict = dict(pool.map(stalling_counting, transcript_list))
    return transcript_dict


def output(transcript_dict):
    outfile = args.output
    with open(outfile, 'w') as out:
        for transcript in seq_dict:
            out.write(transcript + '\n')
            for content in ['seq', 'counts', 'coverage']:
                result = list(map(str, transcript_dict[transcript][content]))
                out.write('\t'.join(result) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='Calculate the stalling counts.')
    parser.add_argument('--threads','-t', required=False, type=int, default=cpu_count(), help='Number of threads to use [default: number of CPUs].')
    parser.add_argument('--fasta','-fa', required=True, type=str, help='FASTA file.')
    parser.add_argument('--bam','-b', required=True, type=str, help='BAM file.')
    parser.add_argument('--output', '-o', required=True, type=str, help='Output file.')
    args = parser.parse_args()

    start_time = time.time()
    seq_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    transcript_dict = extract_transcript_and_counting()
    output(transcript_dict)
    end_time = time.time()
    final_time = round(end_time - start_time, 4)
    print(f'The program finished in {final_time} seconds.')
