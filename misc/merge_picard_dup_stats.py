from glob import glob
import os
import sys
from argparse import ArgumentParser

def parse_command(command):
    """parses inputs"""
    parser = ArgumentParser(description='specific columns from multiple tsv files on a common key column')
    parser.add_argument('-i', '--inputs', action='append', nargs='+', required=True, help='path to folder containing input files')
    #parser.add_argument('-p', '--path', required=True, help='path to folder containing input files')
    parser.add_argument('-o', '--out', required=True, help='output file prefix')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)


def main():
    input_command = parse_command(sys.argv[1:])
    files = input_command.inputs
    prefix = input_command.out
    prefix = 'test_out'

    ## define header
    header = ['SAMPLE','UNPAIRED_READS_EXAMINED','READ_PAIRS_EXAMINED','SECONDARY_OR_SUPPLEMENTARY_RDS',
     'UNMAPPED_READS','UNPAIRED_READ_DUPLICATES', 'READ_PAIR_DUPLICATES', 'READ_PAIR_OPTICAL_DUPLICATES',
     'PERCENT_DUPLICATION', 'ESTIMATED_LIBRARY_SIZE']
    
    # get filenames
    files = glob('/mnt/c/scratch/atac_seq/dedup' + '/*')

    # extract metrics lines and write out to single file
    with open(prefix +".tsv", 'w') as out:
        out.write('\t'.join(header) + '\n')
        for file in files:
            sname = os.path.basename(file).replace('.sorted.duplication_metrics', '')
            print(sname)
            with(open(file)) as f:
                line = f.readlines()[7].split('\t')
                line[0] = sname
                out.write('\t'.join(line) + '\n')

main()
