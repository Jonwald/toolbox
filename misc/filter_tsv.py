import pandas as pd
from glob import glob
import os
import sys
from argparse import ArgumentParser

def parse_command(command):
    """parses inputs"""
    parser = ArgumentParser(description='filter specific rows from a tsv files matching values, eg: for filteting gene expresssion matrices')
    parser.add_argument('-i', '--intsv', required=True, help='path tsv file to filter')
    parser.add_argument('-g', '--genelist', required=True, help='single column file with genes to filter')
    parser.add_argument('-k', '--key', required=True,type=int, help='key column name to match')
    parser.add_argument('-o', '--out', required=True, help='output file prefix')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)


def main():
    # read inputs
    input_command = parse_command(sys.argv[1:])

    # get files and genelist
    tsv = input_command.intsv
    genelist = input_command.genelist
    # key column in tsv 
    key = input_command.key
    # output prefix
    outname = input_command.out

    df = pd.read_csv(tsv, sep='\t')
  
    with open(genelist) as f:
        ids = f.read().splitlines()

    print("filtering")
    output_df = df[df[key].isin(ids)]
    # write output
    output_df.to_csv(outname +".tsv", sep="\t", header=True, index=False, na_rep="NA")

