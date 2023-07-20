import pandas as pd
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

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)

def main():
    input_command = parse_command(sys.argv[1:])
    path = input_command.path
    prefix = input_command.out

    # get files with extension
    stat_files = glob(path + '/*' + ".stats")
    idx_files = glob(path + '/*' + ".idxstats")
        
    # merge glagstat files
    df_list = []
    i=0
    for k in stat_files:
        sample = os.path.basename(k).replace(".stats", '')
        print("reading: " + sample)
        with(open(k)) as f:
            for line in f.readlines():
                if "in total" in line:
                    qc_pass = int(line.split(" ")[0])
                    qc_fail = int(line.split(" ")[2])
                    total = qc_pass + qc_fail
                    print (qc_pass, qc_fail, total)
                if  "primary mapped" in line:
                    primary_mapped = int(line.split(" ")[0])
                    print (primary_mapped)
                map_perc =  round((primary_mapped/total)*100,2)
        df = pd.DataFrame({'sample': sample, 'qc_pass_reads': qc_pass, 'qc_fail_reads': qc_fail, 'total_reads': total, 'primary_mapped_reads': primary_mapped, 'mapped_perc': map_perc}, index=[i])
        i += 1
        df_list.append(df)

    # write out files        
    aln_stats = pd.concat(df_list)
    col_order = ['sample', 'total_reads', 'qc_pass_reads', 'qc_fail_reads', 'primary_mapped_reads', 'mapped_perc']
    aln_stats = aln_stats[col_order]
    aln_stats.to_csv(prefix + "aln_stats.tsv", sep="\t", header=True, index=False, na_rep="NA")
    
    ## merge idxstat files
    df_list = []
    for k in idx_files:
        sample = os.path.basename(k).replace(".idxstats", '')
        print("reading: " + sample)
        df = pd.read_csv(k, sep='\t', names = ['chr', 'length', sample, 'unmapped_reads'], header=None)[['chr', sample]]
        df_list.append(df)
    idx_stats = reduce(lambda x, y: pd.merge(x, y, on = 'chr'), df_list)
    idx_stats.to_csv(prefix + "idx_stats.tsv", sep="\t", header=True, index=False, na_rep="NA")

main()
