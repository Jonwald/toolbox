import csv
import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()
#fname = sys.argv[1]
#taxids = open(fname)

taxids = open("test.csv")

tax_ids2 = csv.reader(taxids)
for taxid in tax_ids2:
	print(ncbi.get_taxid_translator((taxid)))

