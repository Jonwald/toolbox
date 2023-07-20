# meth_perc_report
from glob import glob
import os
import sys
from argparse import ArgumentParser
from functools import reduce


def parse_command(command):
    """parses inputs"""
    parser = ArgumentParser(description='specific columns from multiple tsv files on a common key column')
    parser.add_argument('-p', '--path', required=True, help='path to folder containing input files')
    parser.add_argument('-o', '--out', required=True, help='output file prefix')
    parser.add_argument('-d', '--dp', required=True, type = int, help='depth cutoff to include site')
    parser.add_argument('-m', '--pm', required=True, type = float, help='% methylation cutoff to include site')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)


def summarise_calls(file, cov_cutoff, meth_cutoff, target_chroms):
    with(open(file)) as f:
        counter = {}
        next(f)
        for line in f:
            ll = line.split('\t')
            chr= ll[0]
            if chr not in target_chroms:
                continue
            if chr not in counter.keys():
                counter[chr] = ''
                all_meth = 0
                all_c = 0
            meth = int(ll[4])
            no_meth = int(ll[5])
            cov = meth + no_meth
            if cov < cov_cutoff:
                continue
            if (meth/cov) > meth_cutoff:
                all_meth +=1
            all_c +=1
            counter[chr] = round(all_meth/all_c, 2)
    return(counter)

def main():
    input_command = parse_command(sys.argv[1:])
    path = input_command.path
    prefix = input_command.out
    dp = input_command.dp
    pm = input_command.pm

    target_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY', 'chrM', 'lambda_control', 'pUC19c_control']
    
    perc_meth = summarise_calls(glob(path + '/prefix' + ".bedGraph"), dp, pm, target_chroms)
    for k,v in perc_meth.items():
        if v == '':
            perc_meth[k] = 'NA'
    with open(prefix +".tsv", 'w') as f:
        for k,v in perc_meth.items():
            f.write(k + '\t' + str(v) + '\n')

main()
hdfhbsd
dsfs