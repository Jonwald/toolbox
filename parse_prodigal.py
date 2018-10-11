#!/usr/bin/env python

import gffutils
import argparse

#parse the arguments
parser = argparse.ArgumentParser(description="""Parse Prodigal gff output to produce simple TSV output of contig id to gene ID.""")
parser.add_argument('gff_file', metavar='GFF3', help='Annotation file from Prokka in GFF3 format')
parser.add_argument('--output', default='contig_to_orf.txt', help='Contig ID to ORF ID TSV file (Default: contig_to_orf.txt)')


args = parser.parse_args()

#Input and output files
GFF = args.gff_file
OUT_CDS = open(args.output, 'w')

#load prodigal GFF3 file
db = gffutils.create_db(GFF, ':memory:')

#Print TSV headers
OUT_CDS.write("Contig\tgene_id\tstart\tstop\tstrand\n")


#parse the GFF3 file and write results to output files

for feature in db.all_features():

    start = feature.start - 1
    stop = feature.stop
    seqlen = stop-start

    if (float(start - stop)/float(3)).is_integer() == True:
        partial = str(0)
    else:
        partial = str(1)

    OUT_CDS.write('%s\t%s\t%d\t%d\t%s\n' %(feature.seqid, feature.id, feature.start, feature.stop, feature.strand))
