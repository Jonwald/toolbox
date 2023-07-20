import os
import csv
import pandas as pd
import numpy as np

def read_tsv(file):
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            yield row


def random_gene_names(num_genes):
    """generates random gene names"""
    gene_names = []
    for i in range(num_genes):
        gene_names.append('gene_' + str(i))
    return gene_names



def dummy_expression_matrix(samples, genes):
    """"generates a dummy gene tpm matrix"""
    gene_names=[]
    for i in range(genes):
        gene_names.append('gene_' + str(i))
    sample_names=[]
    for i in range(samples):
        sample_names.append('sample_' + str(i))
    df = pd.DataFrame(np.random.randint(0, 100, size=(genes, samples)), columns=sample_names, index=gene_names)
    return df

def subsample_dataframe(df, fraction):
    """iterate over dataframe columns and downsample to fraction of total counts""" 
    for column in df:
        total_counts = df[column].sum()
        print(total_counts)
        downsampled_counts = round(total_counts * fraction)
        for i in range(downsampled_counts):
            df[column].sample(n=1, replace=True)
            #pd.Index(np.random.choice(df.index, 1, replace=True))
            

    return df

def main():
    test_file = dummy_expression_matrix(10, 100)

    test_file.reset_index(inplace=True)
    #onecolumn_df = pd.melt(test_file, id_vars=['index'])
    # for each colummn
    # sum the counts
    # iterate over columns and downsample to fraction of total counts
    
    ## down
    