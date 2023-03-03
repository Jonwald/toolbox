import random
import re
import sys
from Bio import SeqIO
from argparse import ArgumentParser
from Bio.Seq import Seq
import gzip

def parse_command(command):
    parser = ArgumentParser(description='generates fake pe reads from a fasta input, reads always start from ends of fasta fragment \n'
                                        'readlength and number of fragments and random substiution rate can be controlled \n'
                                        'base quality scores cannot be controlled and are always output as 69 (E) \n'
                                        'R1 always 5 prime , R2 always 3 prime')
    parser.add_argument('-f', '--fasta', required=True, help='input reference fasta')
    parser.add_argument('-l', '--readlen', required=True, help='length of fake paired end reads')
    parser.add_argument('-n', '--numfrags', required=True, help='length of fake paired end reads')
    parser.add_argument('-e', '--err_prob', required=True, help='probabilty of substitution (per base) value between 0.00 and 1')
    parser.add_argument('-o', '--out_prefix', required=True, help='output prefix for output files')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)

input_command = parse_command(sys.argv[1:])


def main():
## inputs ##
    input_file=input_command.fasta
    out_prefix=input_command.out_prefix
    err_prob=input_command.err_prob
    readlen=input_command.readlen
    numfrags=input_command.numfrags

    lib = "ATGC"
    ref = SeqIO.read(input_file,'fasta')


    output_read1 = out_prefix + "_" + str(int(err_prob * 100)) + "_R1.fastq.gz"
    output_read2 = out_prefix + "_" + str(int(err_prob * 100)) + "_R2.fastq.gz"
    output_handle = gzip.open(output_read1 , "wt")
    output_handle2 = gzip.open(output_read2, "wt")

    for i in range(numfrags):
    # loop over x times
    # for each char in fasta, if rand prob > error_prob, then add snp
        err_ref = ''.join([x if random.random() >= err_prob else random.choice(re.sub(x, "", lib)) for x in list(ref.seq)])
        r1_id = ref.id + "_5p_" + str(int(err_prob * 100)) + "_" + str(readlen)
        r2_id = ref.id + "_3p_" + str(int(err_prob * 100)) + "_" + str(readlen)
        r1_seq = err_ref[:readlen]                                 # generate R1
        r2_seq = str(Seq(err_ref).reverse_complement()[:readlen])  # generate R2
        r1_qual = "E" * readlen
        r2_qual = "E" * readlen
        output_handle.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_qual))
        output_handle2.write("@%s\n%s\n+\n%s\n"% (r2_id, r2_seq, r2_qual))

    output_handle.close()
    output_handle2.close()


main()
