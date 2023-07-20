## requires python 3+

## requires following python libraries: ete3 and collections

## takes a list of ncbi taxids, finds the entire lineage and prints the results

## Note; this script vurrently outputs the desired ranks but columns are NOT in the specified order,

## final tsv file will need to be manually reordered (or using awk / sed) until a fix can be found for this behaviour

import csv
from ete3 import NCBITaxa
from collections import OrderedDict
import sys

ncbi = NCBITaxa()
fname = sys.argv[1]
taxids = open(fname)

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = OrderedDict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}


if __name__ == '__main__':
    desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    results = list()
    for taxid in taxids:
        results.append(list())
        results[-1].append(str(taxid))
        ranks = get_desired_ranks(taxid, desired_ranks)
        for key, rank in ranks.items():
            if rank != '<not present>':
                results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
            else:
                results[-1].append(rank)


#generate the header
header = ['Original_query_taxid']
header.extend(desired_ranks)
print('\t'.join(header))

#print the results
for result in results:
	print('\t'.join(result))
