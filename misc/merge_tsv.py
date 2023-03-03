import pandas as pd
from glob import glob
import os
import sys
from argparse import ArgumentParser

def parse_command(command):
    """parses inputs"""
    parser = ArgumentParser(description='specific columns from multiple tsv files on a common key column')
    parser.add_argument('-p', '--path', required=True, help='path to folder containing files to merge')
    parser.add_argument('-e', '--ext', required=True, help='extension of files to merge eg: ReadsPerGene.out.tab')
    parser.add_argument('-k', '--key', required=True, help='key column to merge on')
    parser.add_argument('-c', '--col', required=True, help='value column to keep')
    parser.add_argument('-g', '--genelist', required=True, help='single column file with genes to filter')
    parser.add_argument('-o', '--out', required=True, help='output file prefix')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)

def main():
    input_command = parse_command(sys.argv[1:])

    # get files with extension
    path = input_command.path
    ext = input_command.ext

    # key to match on
    key = input_command.key
    # column to merge
    col = input_command.col
    outname = input_command.out
    #genelist
    genelist = input_command.genelist
    # get file names in list
    tsv_files = glob(path + '/*' + ext)

    if len(tsv_files) < 2:
        print("error, requires at least two files to merge")
        exit

    with open(genelist) as f:
        ids = f.read().splitlines()

    # loop over files and merge on key
    df_list = []
    i=0
    for k in tsv_files:
        print("reading: " + k)
        tsv = pd.read_csv(k, sep='\t')[[key,col]]
        tsv = tsv[tsv[key].isin(ids)]
        name = os.path.basename(k).replace(ext, '')
        tsv.rename(columns={key: "GENE", col: name}, inplace=True)
        tsv.set_index("GENE", inplace=True)
        tsv = tsv.groupby(tsv.index).agg({name : sum})
        df_list.append(tsv)
        i += 1

    out_df = pd.concat(df_list, axis=1)

    out_df.to_csv(outname +".tsv", sep="\t", header=True, index=True, na_rep="NA")

main()
