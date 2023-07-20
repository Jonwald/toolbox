from os import chdir, getcwd
from pysam import VariantFile
import pysam

def read_bed(bed):
    """reads in .bed file, outputs bed_list, list of dicts
        dict struct {Keys, Start, End, Desc}"""
    #define lines - individual line, bed_list - list of dicts, keys = dict keys
    lines = []
    bed_list = []
    keys = ['Chrom', 'Start', 'End', 'Desc']
    with open(bed, 'r') as file:
        #read in lines from file to lines
        for line in file:
            lines.append(line.strip().split())
        #combine ea row to key, generate dict, append to bed_list
        for row in lines:
            zip_row = zip(keys, row)
            dict_row = dict(zip_row)
            bed_list.append(dict_row)

    #convert bed_list "Start" and "End" to int() for comparison
    for region in bed_list:
        region["Start"] = int(region["Start"])
        region["End"] = int(region["End"])

    return bed_list


in_file = 'test.vcf'
with open(in_file, 'r') as vcf_in:
    vf = VariantFile(in_file, 'r')
    vf.fetch()
            if rec.info["DP"] >= 10:
                print(rec.info)

def get_rr(path, depth):
    vcf_in = VariantFile(path)
    for rec in vcf_in.fetch():
        if rec.info["DP"] >= depth:




    
        



def vcf_variant_parse2(in_file, bed_list, out_path, sample_id, info_id):
    """Takes an input vcf file, filters out variants, defines list of variants for downstream qc comparison
            in_file -- filepath to vcf file being analyzed
    """
    #read in_file; define filter_file output to filter to
    with open(in_file, 'r') as vcf_file:
        vcf_orig = vcf.Reader(vcf_file)
        vcf_orig.metadata["Sample ID"] = sample_id
        vcf_orig.metadata["Info ID"] = info_id
        out_vcf = os.path.join(out_path, sample_id + "_" + timestamp + "_filtered.vcf")
        #iterate over in_file to filter variants to filter_file
        list_variant = []
        samp_id = vcf_orig.samples[0]
        with open(out_vcf, 'w') as filtered_vcf:
            # write header
            write_vcf = vcf.Writer(filtered_vcf, vcf_orig)
            for record in vcf_orig:
                vcf_orig.metadata["Sample ID"] = sample_id
                vcf_orig.metadata["Info ID"] = sample_id
                #define call for VF selection
                call = record.genotype(samp_id)
                mask = False
                #applied filters: variant call detected, FILTER == 'PASS', VF >= 3.0%
                if record.ALT != [None] and record.FILTER == [] and call.data.VF > 0.03:
                    #apply demask, only include variant if in bed region
                    for region in bed_list:
                        if record.CHROM == region["Chrom"] and record.POS in range(region["Start"], region["End"]):
                            mask = True
                    for sample in record.samples:
                        var_freq = sample['VF']
                    if 'ANT' in list(record.INFO.keys()):
                        var_ann = record.INFO['ANT']
                    if 'ANT' not in list(record.INFO.keys()):
                        var_ann = ''

                if record.ALT != [None] and record.FILTER == [] and call.data.VF >= 0.03 and mask == False:
                    variant = [record.CHROM, record.POS, record.REF, str(record.ALT[0]), var_freq, var_ann]
                    list_variant.append(variant)
                    write_vcf.write_record(record)
    return list_variant




def get_rr_cov(vcf_file, rr_file):
    with open(vcf_file, 'r') as vcf_file:
        vcf_orig = vcf.Reader(vcf_file)
            #iterate over in_file to filter variants to filter_file
        d = defaultdict(int)
        for record in vcf_orig:
            coord = record.CHROM + "_" + str(record.POS)
            d[coord] += record.INFO['DP']
        d = dict((k, v) for k, v in d.items() if v >= 259)
        covered_pos = [x for x in dict.keys(d)]
        with open(rr_file) as r:
            rr_reader = csv.reader(r, delimiter="\t")
            reportable_range = [tuple(row) for row in rr_reader]
            out = [item for t in reportable_range for item in t]

            n_covered = len(list(set(covered_pos) & set(out)))
            rr_covered = (n_covered / (len(set(out))-1) * 100)
            if rr_covered >= 95:
                result = True
            else:
                result = False
    return result