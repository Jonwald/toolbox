#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 13:00:51 2020

@author: jbyoung
"""

import sys
import csv
from argparse import ArgumentParser
import numpy
from cigar import Cigar


def parse_command(command):
    """Sets up and parses input arguemnts"""
    parser = ArgumentParser(description='Counts number of unique start sites for each fusion \n'
                                        'in star fusion results, takes a list of star fusion \n'
                                        'results folders as input, usage:\n'
                                        'get_uniq_starts_final.py -f *star_fusion_results -o output_name.txt')
    parser.add_argument('-f', '--fusionresults', nargs='+', required=True, help='folder(s) containing star-fusion results')
    parser.add_argument('-o', '--outname', required=True, help='output for combined results file')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return args

def chomp_cigar_list(table, gene):
    """Reads cigar string for supporting read pairs and calculates start/end points"""
    coord = []
    if gene == "A":
        aln_start = 0
        splice_junc = 1
        strand_ind = 2
        cigar_ind = 3
    if gene == "B":
        aln_start = 4
        splice_junc = 5
        strand_ind = 6
        cigar_ind = 7
    for row in table:
        frag_type = row[9]
        #print("fragtype = " + frag_type)
        c_string = Cigar(row[cigar_ind].replace("-", ""))
        c_list = list(c_string.items())
        c_ar = numpy.array(c_list)
        p_ind = numpy.where(c_ar == 'p')

        if p_ind[0].size > 0:
            c_ar = c_ar[:int(p_ind[0])]

        clipped = numpy.where(c_ar == 'S')
        if clipped[0].size > 0:
            c_ar = numpy.delete(c_ar, clipped[0], 0)
        deletions = numpy.where(c_ar == 'D')
        if deletions[0].size > 0:
            c_ar = numpy.delete(c_ar, deletions[0], 0)

        mapped_len = numpy.sum(c_ar[:, 0].astype(numpy.int))

        if gene == "A" and row[strand_ind] == '+':
            if frag_type == "spanningfrag":
                coord.append(int(row[aln_start]))
            else:
                coord.append(int(row[splice_junc]) - int(mapped_len))

        if gene == "A" and row[strand_ind] == '-':
            if frag_type == "spanningfrag":
                coord.append(int(row[aln_start]) + int(mapped_len))
            else:
                coord.append(int(row[splice_junc]) + int(mapped_len))

        if gene == "B" and row[strand_ind] == '+':
            if frag_type == "spanningfrag":
                coord.append(int(row[aln_start]) + int(mapped_len))
            else:
                coord.append(int(row[splice_junc]) + int(mapped_len))

        if gene == "B" and row[strand_ind] == '-':
            if frag_type == "spanningfrag":
                coord.append(int(row[aln_start]))
            else:
                coord.append(int(row[splice_junc]) - int(mapped_len))
    return coord


def decomment(csvfile):
    """strip comment lines from csv"""
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw:
            yield raw

def cluster_ends(coords, dist):
    """groups read ends < dist bp apart as same primer"""
    ref = 0
    cnt = 0
    out = []
    for pos in coords:
        if abs(pos[1]-ref) > dist:
            cnt += 1
            ref = pos[1]
        out.append((str(pos[0]), 'primer%d'%cnt))
    return out


def main():
    """main function for counting unique starts from star fusion supporting reads"""
    input_command = parse_command(sys.argv[1:])
    csv.field_size_limit(sys.maxsize)

    with open(input_command.outname, "w+") as out:
    #with open("test_out", "w+") as out:
        writer = csv.writer(out, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['Sample', 'Fusion', 'JunctionReads', 'SpanningFragments', 'UniqueStarts'])
        for path in input_command.fusionresults:
        #for path in ["S1568-00245-CP1Lib1_S23_Resampled_uniq.R1.fusions"]:
            fusions_file = str(path) + "/star-fusion.fusion_predictions.tsv"
            coding_eff_file = str(path) + "/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
            chimeric_reads_file = str(path) + "/Chimeric.out.junction"
            print("starting: " + path)
        # need fields
        # 0: fusion name,
        # 1: junction reads
        # 2: spanning frags
        # 8: junction readnames
        # 9: spanningfragnames

        # dictionary:
        # readname: s_start, a_junction, a_strand, a cigar, b_start, b_junction, b_strand, b_cigar

            junc_dict = {}
            with open(chimeric_reads_file, newline='') as chim_in:
                next(chim_in)
                reader = csv.reader(decomment(chim_in), delimiter='\t')
                for row in reader:
                    read = row[9]
                    stats = [row[10], row[1], row[2], row[11], row[12], row[4], row[5], row[13], row[9]]
                    junc_dict[read] = stats

                with open(fusions_file, newline='') as fus_in, open(coding_eff_file, newline='') as eff_in:
                    next(fus_in)
                    next(eff_in)
                    results = [x for x in csv.reader(fus_in, delimiter='\t')]
                    eff = [x for x in csv.reader(eff_in, delimiter='\t')]
                for k in range(0, (len(results))):
                    line=results[k]
                #for line in results:
                    splitfrags = line[8].split(',')
                    spanningfrags = line[9].split(',')
                    reads = (splitfrags + spanningfrags)
                    reads = list(filter(lambda a: a != '.', reads))
                    read_pos = list((junc_dict[k]) for k in tuple(reads))
                    for i in read_pos:
                        if i[8] in spanningfrags:
                            i.append("spanningfrag")
                        else:
                            i.append("splitfrag")
                    a_coords = chomp_cigar_list(read_pos, "A")
                    b_coords = chomp_cigar_list(read_pos, "B")
                    res = list(zip(a_coords, b_coords))
                    comb_list2 = sorted(res, key=lambda x: x[1])
                    out = cluster_ends(comb_list2, 10)
                    clus = len(set(out))
                    coding_effect = eff[k][19]
                    writer.writerow([path.replace("_uniq.R1.fusions", ""), line[0], line[1], line[2], clus, coding_effect])
                    print("finished: " + path + " " + line[0])

main()

