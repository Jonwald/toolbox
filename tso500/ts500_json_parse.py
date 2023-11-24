## Parse TSO JSON file
import sys
import os
from argparse import ArgumentParser
import collections
from itertools import combinations
import json
import csv
import gzip

"""
old headings
"Gene, Chromosome, Genomic Position, Reference Call, 
Alternative Call, Allele Frequency, Depth, 
P-Dot Notation, C-Dot Notation, Consequence(s), Affected Exon(s)"

new headings
Gene, Chromosome, Genomic Position, Reference Call, 
Alternative Call, Allele Frequency, Depth, 
P-Dot Notation, C-Dot Notation, Consequence(s), Affected Exon(s),
Combined allAf, Somatic/Germline Call
"""


"""
filtering logic

 Somatic if
Annotated as “Missense”, “frameshift”, “nonsense”, “stop_gained” or “start_lost”
AND  
Combined gnomad and gnomadExome “allAf” < 0.01
AND
[
    COSMIC count >= 20
    OR
    [
        Variant frequency<0.90 
        AND 
        Not classed as “inherited”, “paternal”, “maternal”, “biparental”, “uniparental” by ClinVar 
    ]
]

Germline if
None of the criteria above are met. """


## functions

def parse_command(command):
    """Parses command line arguments"""
    parser = ArgumentParser(description='modify TSO500 combinedvariantoutput file to add allAF, somatic/germline call and functional significance \
                                         columns and write out to a new file.')
    parser.add_argument('-i', '--input', required=True, help='input folder path, containing "annotated json and combinedvariantoutput files')
    parser.add_argument('-o', '--output', required=True, help='output folder path, to write the updated combinedvariantoutput files to')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)

    
def parse_tmb(tmb_file):
    """
    Read the TMB file and return a dict with the max cosmic counts using chr_pos_ref_alt as key
    """
    with open(tmb_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        tmb_dict = {}
        for row in reader:
            chrom = row['Chromosome']
            pos = row['Position']
            ref = row['RefCall']
            alt = row['AltCall']
            index = '_'.join([chrom, pos, ref, alt])
            tmb_dict[index] = int(row['MaxCosmicCount'])
    return tmb_dict
    
def read_json(j_file, tmb_dict):
    """
    Read the TSO500 JSON file and return a dict with the gnomad counts 
    and clinvar annotations using chr_pos_ref_alt as key
    """
    with gzip.open(j_file, 'rt', encoding='UTF-8') as f:
        data = json.load(f)

    # create dict with chr_pos_ref_alt as key and gnomad counts + clinvar annotations as values
    af_dict = {}       
    for i in (x for x in data['positions'] if x['filters'][0] == 'PASS'):
        if "altAlleles" not in i.keys():
            continue
        chrom = i['chromosome']
        pos = i['position']
        ref = i['refAllele']
        alt = i['altAlleles'][0]
        index = '_'.join([chrom, str(pos), ref, alt])
        
        # get af counts and calc allAF
        variant = i['variants'][0]
        gn_ac = variant.get('gnomad', {}).get('allAc', 0)
        gn_an = variant.get('gnomad', {}).get('allAn', 0)
        ex_ac = variant.get('gnomadExome', {}).get('allAc', 0)
        ex_an = variant.get('gnomadExome', {}).get('allAn', 0)

        try:
            all_af = round((gn_ac + ex_ac)/(gn_an + ex_an), 6)
        except ZeroDivisionError:
            all_af = 0
        
        # get max cosmic count
        # use lambda function to get highest sample count in list of dicts
        cos_count = tmb_dict[index]
              
        # get clinvar annotations if they exist
        clinvar = [item for sublist in [x['alleleOrigins'] for x in variant.get('clinvar', [])] for item in sublist] or ["NA"]
        
        # populate dict
        af_dict[index] = [all_af, cos_count, clinvar]
        
    return af_dict

def get_sv_start(v_file):
    """
    Find the line wherethe small variants table starts in the combined variant output file.
    """
    with open(v_file, 'r') as f:
        # find line number which starts with [Small Variants]
        for num, line in enumerate(f):
            if line.startswith('[Small Variants]'):
                table_start = num + 1
                break
    return table_start

def is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf):
    """
    Determine if a variant is somatic based on the following criteria:
    AF < 0.01 and (cosmic count >= 20 or (AF < 0.9 and no clinvar germline annotations))
    """
    if any(x == y for x in fun_ann for y in cons_filter) and \
        af_dict[index][0] < 0.01 and \
            (af_dict[index][1] >= 20 or \
                (vf < 0.9 and \
                    not any(x == y for x in af_dict[index][2] for y in clin_ann))):
        return True
    else:
        return False

def get_funsig(af_dict, index):
    """
    return functional significance of a variant based on the following criteria:
    'Known': Cosmic count >= 20 & Combined gnomad and gnomadExome “allAf” < 0.001
    'Likely': (truncating mutation or Cosmic count >= 5 & < 20) & ( Combined gnomad and gnomadExome “allAf” < 0.001)
    'Unknown': All remaining calls
    TODO change trunc mutation to anything ending in TER in the p dot notation
    not to be used for now
    """
    tunc_mutations = ['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'stop_lost', 'start_lost']
    
    all_af = af_dict[index][0]
    cosmic_count = af_dict[index][1]
    mutations = af_dict[index][2]
    
    is_known = cosmic_count >= 20 and all_af < 0.001
    is_likely = (any(mutation in mutations for mutation in tunc_mutations) or 5 <= cosmic_count < 20) and all_af < 0.001

    if is_known:
        return 'Known'
    elif is_likely:
        return 'Likely'
    else:
        return 'Unknown'

def copy_first_x_lines(input_file, output_file, x):
    """
    copy the top part of the input_file up to linenum x
    and write out to output_file
    """
    with open(input_file, 'r') as in_file, open(output_file, 'w', newline='') as out_file:
        reader = csv.reader(in_file, delimiter='\t')
        writer = csv.writer(out_file, delimiter='\t')

        for i, row in enumerate(reader):
            if i == x:
                break
            writer.writerow(row)

def get_file_pairs(directory):
    """
    Find all pairs of files ending with .json.gz or CombinedVariantOutput.tsv
    """
    # Get all file names in the directory
    filenames = os.listdir(directory)
    relevant_files = [filename for filename in filenames if filename.endswith('.json.gz') or filename.endswith('CombinedVariantOutput.tsv') or filename.endswith('TMB_Trace.tsv')]
    # Get the prefixes (everything before the first underscore)
    prefixes = [filename.split('_')[0] for filename in relevant_files]
    prefixes = [prefix.replace('-D1', '') for prefix in prefixes]

    # Group filenames by prefix
    prefix_dict = collections.defaultdict(list)
    for prefix, filename in zip(prefixes, relevant_files):
        if filename.endswith('.json.gz') or filename.endswith('CombinedVariantOutput.tsv') or filename.endswith('TMB_Trace.tsv'):
            prefix_dict[prefix].append(filename)
    # convert dict entires to sets
    for k,v in prefix_dict.items():
        prefix_dict[k] = set(v)

    return prefix_dict


def main():
    ## consequences / annotations to filter on   
    
    ## 7th sept
    """cons_filter = ['missense_variant', 'frameshift', 'nonsense', 'stop_gained', 'start_lost']"""
    
    ## 16th oct
    """cons_filter = ['missense_variant', 'frameshift_variant', 'stop_gained', 'start_lost']"""

    ## corrected
    """cons_filter = ['missense_variant', 'frameshift_variant','splice_acceptor_variant', 
                            'splice_donor_variant', 'inframe_deletion', 'inframe_insertion',
                            'synonymous_variant', 'stop_gained', 'stop_lost', 'start_lost',
                            'upstream_gene_variant']"""
    
    ## 20th Nov
    cons_filter = ['missense_variant', 'frameshift_variant', 'splice_acceptor_variant',
                     'splice_donor_variant', 'inframe_deletion','inframe_insertion',
                     'stop_gained', 'stop_lost', 'start_lost', '3_prime_UTR_variant',
                     '5_prime_UTR_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                     'transcript_ablation', 'transcript_amplification', 'feature_elongation',
                      'feature_truncation', 'protein_altering_variant']

    ## clinvar germline annotations
    clin_ann = ['inherited', 'paternal', 'maternal', 'biparental', 'uniparental']

    ## read an input directory and find required files
    parser = parse_command(sys.argv[1:])
    datadir = parser.input
    outdir = parser.output
    file_pairs = get_file_pairs(datadir)

    for prefix, pair in file_pairs.items():
        j_file = os.path.join(datadir, [x for x in pair if x.endswith('.json.gz')][0])
        v_file = os.path.join(datadir, [x for x in pair if x.endswith('CombinedVariantOutput.tsv')][0])
        tmb_file = os.path.join(datadir, [x for x in pair if x.endswith('TMB_Trace.tsv')][0])
        print("processing: " + prefix)

        # make a dict with the gnomad counts using chr_pos_ref_alt as key
        print("reading: " + tmb_file)
        tmb_dict = parse_tmb(tmb_file)

        print("reading: " + j_file)
        af_dict = read_json(j_file, tmb_dict)
    
        # read combined_variant_out
        print("reading: " + v_file)
       
        table_start = get_sv_start(v_file)

        ## set up output file
        output_file = os.path.join(outdir, prefix + '_CombinedVariantOutputwithGermlinecall.tsv') 
        copy_first_x_lines(v_file, output_file, table_start)

        with open(v_file, 'r') as f, open(output_file, 'a', newline='') as out_file:
            reader = csv.DictReader(f.readlines()[table_start:], delimiter='\t')
            outfields = reader.fieldnames + ['allAF', 'Somatic/Germline Call']
            writer = csv.DictWriter(out_file, fieldnames=outfields, delimiter='\t')
            
            ## write header to outfile
            writer.writeheader()
            print("applying logic and writing output: " + prefix + '_CombinedVariantOutputwithGermlinecall.tsv')
            somatic_counts = 0
            germline_counts = 0
            for row in reader:
                if not row['Depth']:  
                    continue
                index = '_'.join([row['Chromosome'], str(row['Genomic Position']), row['Reference Call'], row['Alternative Call']])
                fun_ann = row['Consequence(s)'].split(':')

                ## get variant frequency 
                vf = float(row['Allele Frequency'])

                # apply logic to get classification
                if is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf):
                    sg_call = "Somatic"
                    somatic_counts += 1
                else:
                    sg_call = "Germline"
                    germline_counts += 1
                                    
                # funsig = get_funsig(af_dict, index)
                
                ## add new columns to row
                row['allAF'] = af_dict[index][0]
                row['Somatic/Germline Call'] = sg_call
                # row['Functional significance'] = funsig
                ## write row to out_file
                writer.writerow(row)
            print("total somatic variants: " + str(somatic_counts))
            print("total germline variants: " + str(germline_counts))
            
if __name__ == '__main__':
    main()
