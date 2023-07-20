import pysam
from sys import argv
import csv

bam_file = argv[1]
bed_file = argv[2]
out_prefix = argv[3]

#bam_file = "/home/jbyoung/toolbox/toolbox/bam/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam"
#bed_file = "/home/jbyoung/toolbox/toolbox/bam/all_chr_windows"

bam = pysam.AlignmentFile(bam_file, 'rb')

with open(bed_file) as bf:
    row_out=[]
    with open(out_prefix + ".csv", 'w') as file:
        for line in bf:
            line_parts = line.strip().split()
            chr = line_parts[0]
            start = int(line_parts[1])
            end = int(line_parts[2])
            read_data = bam.fetch(chr, start, end)
            total_reads = 0
            total_bases = 0
            gc_bases = 0
            for read in read_data:
                total_reads += 1
                seq = read.query_sequence
                total_bases += len(seq)
                gc_bases += len([x for x in seq if x == 'C' or x == 'G'])
            if total_bases == 0:
                gc_percent = 'No Reads'
            else:
                gc_percent = float(gc_bases)/total_bases * 100
            row_out.append([chr, total_reads, gc_percent])
        write = csv.writer(file)
        write.writerows(row_out)
        #print '{0}\t{1}'.format(line.strip(), gc_percent)