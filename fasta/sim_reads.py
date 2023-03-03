from distutils.log import error
from ftplib import error_proto
import random
import re
import os
import sys
from xml.dom.expatbuilder import FragmentBuilderNS
from Bio import SeqIO
from argparse import ArgumentParser
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import gzip

def parse_command(command):
    parser = ArgumentParser(description='models PE reads over given target regions for a given refernec genome\n'
                                        'no of fragments, average length and stdev of fragments can be coneolled \n'
                                        'per base random substiution rate for each read can be controlled \n'
                                        'base quality scores cannot be controlled and are always output as 69 (E) \n'
                                        'assumes random fragmentation over the target region \n'
                                        'target region can toptionally be written to fasta \n'
                                        'R1 always fw strand, R2 always reverse stand, and generated from the ends of the fragments')
    parser.add_argument('-f', '--fasta', required=True, help='input reference fasta')
    parser.add_argument('-t', '--targets', required=True, help='input target regions in tsv format')
    parser.add_argument('-s', '--snps', required=True, help='input snps to model in tsv format')
    parser.add_argument('-l', '--readlen', required=True, type=int, help='length of modelled paired end reads')
    parser.add_argument('-n', '--numreads', required=True, type=int, help='number of read pairs to generate')
    parser.add_argument('-e', '--err_prob', required=True, type=float, help='probabilty of substitution (per base) value between 0.00 and 1')
    parser.add_argument('-fl', '--fraglen', required=True, type=int, help='average fragment length')
    parser.add_argument('-sd', '--stdev', required=True, type=int, help='standard devation of fragemnt length')
    parser.add_argument('-o', '--out_prefix', required=True, help='output prefix for output files')
#    parser.add_argument('-wt', type=bool, action='store_true', help='optional, boolean, write out target fasta sequences')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)

input_command = parse_command(sys.argv[1:])


def snp_replace(s, position, character):
    return s[:position] + character + s[position+1:]

def model_snps(record_dict, targets, snps):
    haplo_list=[]
    #target = targets.values[1]
    for target in targets.values:
        seq_1 = record_dict[str(target[0])].seq[target[1]:target[2]]
        seq_2 = record_dict[str(target[0])].seq[target[1]:target[2]]
        variants = snps[snps.amplicon == target[3]] #[target[3]]
        variants['adjusted_coord'] = variants.coord - target[1] - 1
        for var in variants.values:
            genotype= var[5]
            coord = var[7]
            ref = var[3]
            alt = var[4]
            if genotype == 'hom_ref':
                seq_1 = snp_replace(seq_1, coord, ref)
                seq_2 = snp_replace(seq_2, coord, ref)
            if genotype == 'het':
                seq_1 = snp_replace(seq_1, coord, ref)
                seq_2 = snp_replace(seq_2, coord, alt)
            if genotype == 'hom_alt':
                seq_1 = snp_replace(seq_1, coord, alt)
                seq_2 = snp_replace(seq_2, coord, alt)
        haplo_list.append(seq_1)
        haplo_list.append(seq_2)
    return haplo_list

def generate_frags(target_seqs, read_target, model, av_frag_len, sd):
    frag_seq = []
    frag_target = read_target / 2
    frags_per_target = round(frag_target / len(target_seqs))
    
    for seq in target_seqs:
        if model == "gauss":
            frag_lens = np.round(np.random.normal(loc=av_frag_len, scale=sd, size=frags_per_target))
        if model == "random":
            frag_lens = [random.randrange(av_frag_len - sd, av_frag_len + sd, 1) for i in range(frags_per_target)]
        for i in frag_lens:
            start = random.randint(0, len(seq)-int(i))
            frag_seq.append(seq[start:start+int(i)])
    return frag_seq

def write_reads(frags, readlen, err_prob, out_prefix):
        output_read1 = out_prefix + "_" + str(int(err_prob * 100)) + "_R1.fastq.gz"
        output_read2 = out_prefix + "_" + str(int(err_prob * 100)) + "_R2.fastq.gz"
        output_handle = gzip.open(output_read1, "wt")
        output_handle2 = gzip.open(output_read2, "wt")
        i=1
        for frag in frags:
            # get reads
            r1_seq=str(frag[:readlen])
            r2_seq=str(frag[-readlen:].reverse_complement())

            ## model errors
            if err_prob > 0:
                r1_seq = ''.join([x if random.random() >= err_prob else random.choice(re.sub(x, "", "ATGC")) for x in r1_seq])
                r2_seq = ''.join([x if random.random() >= err_prob else random.choice(re.sub(x, "", "ATGC")) for x in r2_seq])
            
            # model qual string -- tba

            # generate header and qual string

            r1_id = "read_" + str(i) + "_1"
            r2_id = "read_" + str(i) + "_2"
            r1_qual= "E" * readlen
            r2_qual= "E" * readlen
            
            output_handle.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_qual))
            output_handle2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_qual)) 
            i+=1
        output_handle.close()
        output_handle2.close()

def main():
    # define inputs
    input_command = parse_command(sys.argv[1:])
    reference_genome=input_command.fasta
    numreads=input_command.numreads
    targets=input_command.targets
    snps=input_command.snps
    out_prefix=input_command.out_prefix
    err_prob=input_command.err_prob
    readlen=input_command.readlen
    fraglen=input_command.fraglen
    sd=input_command.stdev
    # write_targets=input_command.write_targets
    
    # load genome
    print("reading genome")
    record_dict=SeqIO.to_dict(SeqIO.parse(reference_genome, "fasta"))
    
    # load targets and snps of interest
    print("reading targets and snps")
    targets=pd.read_table(targets, delimiter="\t")
    snps=pd.read_table(snps, delimiter="\t")
   
    # model snps
    print("generating haplotypes")
    haplo_targets = model_snps(record_dict, targets, snps)
    write_targets = True
    if write_targets:
        print("writing targets")
        i=1
        output_handle = open("targets_out.fasta", "wt")
        for seq in haplo_targets:
            output_handle.write(">%s\n%s\n" % ("ref_" + str(i), str(seq)))
            i+=1
        output_handle.close()
    
    # generate fragments
    print("generating fragments")
    frags = generate_frags(haplo_targets, numreads, "random", fraglen, sd)

    # write out reads
    print("writing reads")
    write_reads(frags, readlen, err_prob, out_prefix)

main()
