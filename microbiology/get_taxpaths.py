import csv
import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()
#fname = sys.argv[1]


#taxids = open(fname)
taxids = open("test.csv")

for taxid in taxids:
	lineage=ncbi.get_lineage(taxid)
	names = ncbi.get_taxid_translator(lineage)
	print([taxid],[names[taxid] for taxid in lineage])
	

