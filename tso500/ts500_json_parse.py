## Parse TSO JSON file
import sys
import os
from argparse import ArgumentParser
import collections
from itertools import combinations
import json
import csv

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

def read_json(j_file):
    """
    Read the TSO500 JSON file and return a dict with the gnomad counts 
    and clinvar annotations using chr_pos_ref_alt as key
    """
    with open(j_file, 'r') as f:
        data = json.load(f)

    # create dict with chr_pos_ref_alt as key and gnomad counts + clinvar annotations as values
    af_dict = {}       
    for i in (x for x in data['positions'] if x['filters'][0] == 'PASS'):
        chr = i['chromosome']
        pos = i['position']
        ref = i['refAllele']
        alt = i['altAlleles'][0]
        index = '_'.join([chr, str(pos), ref, alt])

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
        try:
            cos_count = max(variant.get('cosmic', []), key=lambda x: x['sampleCount'])['sampleCount']
        except:
            cos_count = 0
              
        # get clinvar annotations if they exist
        clinvar = [item for sublist in [x['alleleOrigins'] for x in variant.get('clinvar', [])] for item in sublist] or ["NA"]
        
        # populate dict
        af_dict[index] = [all_af, cos_count, clinvar]
        
    return af_dict

def get_sv_start(v_file):
    """
    Find the start of the small variants table in the combined variant output file.
    """
    with open(v_file, 'r') as f:
        # find line number which starts with [Small Variants]
        for num, line in enumerate(f):
            if line.startswith('[Small Variants]'):
                table_start = num + 1
                break
    return table_start

def is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann):
    """
    Determine if a variant is somatic based on the following criteria:
    AF < 0.01 and (cosmic count >= 20 or (AF < 0.9 and no clinvar germline annotations))
    """
    if any(x == y for x in fun_ann for y in cons_filter) and \
        af_dict[index][0] < 0.01 and \
            (af_dict[index][1] >= 20 or \
                (af_dict[index][0] < 0.9 and \
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
    Find all pairs of files ending with .json or CombinedVariantOutput.tsv
    """
    # Get all file names in the directory
    filenames = os.listdir(directory)

    # Get the prefixes (everything before the first underscore)
    prefixes = [filename.split('_')[0] for filename in filenames if filename.endswith('.json') or filename.endswith('CombinedVariantOutput.tsv')]

    # Group filenames by prefix
    prefix_dict = collections.defaultdict(list)
    for prefix, filename in zip(prefixes, filenames):
        prefix_dict[prefix].append(filename)

    # Find all pairs of files ending with .json or CombinedVariantOutput.tsv
    file_pairs = [pair for files in prefix_dict.values() for pair in combinations(files, 2)]

    return file_pairs


def main():
    ## consequences / annotations to filter on   
    
    ## 7th sept
    """cons_filter = ['missense_variant', 'frameshift', 'nonsense', 'stop_gained', 'start_lost']"""
    
    ## 16th oct
    cons_filter = ['missense_variant', 'frameshift_variant', 'stop_gained', 'start_lost']

    ## corrected
    """cons_filter = ['missense_variant', 'frameshift_variant','splice_acceptor_variant', 
                            'splice_donor_variant', 'inframe_deletion', 'inframe_insertion',
                            'synonymous_variant', 'stop_gained', 'stop_lost', 'start_lost',
                            'upstream_gene_variant']"""

    ## clinvar germline annotations
    clin_ann = ['germline', 'inherited', 'paternal', 'maternal', 'biparental', 'uniparental']

    ## read an input directory and find required files
    parser = parse_command(sys.argv[1:])
    datadir = parser.input
    outdir = parser.output
    file_pairs = get_file_pairs(datadir)

    for pair in file_pairs:
        j_file = os.path.join(datadir, pair[0])
        v_file = os.path.join(datadir, pair[1])
        
        # make a dict with the gnomad counts + clinvar annotations using chr_pos_ref_alt as key
        print("reading: " + j_file)
        af_dict = read_json(j_file)

        # read combined_variant_out
        print("reading: " + v_file)
       
        table_start = get_sv_start(v_file)

        ## set up output file
        output_file = os.path.join(outdir, v_file.split('.')[0] + 'withGermlinecall.tsv') 
        copy_first_x_lines(v_file, output_file, table_start)

        with open(v_file, 'r') as f, open(output_file, 'a', newline='') as out_file:
            reader = csv.DictReader(f.readlines()[table_start:], delimiter='\t')
            outfields = reader.fieldnames + ['allAF', 'Somatic/Germline Call', 'Functional significance']
            writer = csv.DictWriter(out_file, fieldnames=outfields, delimiter='\t')
            
            ## write header to outfile
            writer.writeheader()
            print("applying logic and writing output: " + v_file.split('.')[0] + 'withGermlinecall.tsv')
            for row in reader:
                if not row['Depth']:  
                    continue
                index = '_'.join([row['Chromosome'], str(row['Genomic Position']), row['Reference Call'], row['Alternative Call']])
                fun_ann = row['Consequence(s)'].split(':')
                
                # apply logic to get classification
                sg_call = "NA"
                
                if is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann):
                    sg_call = "Somatic"
                    print(sg_call)
                else:
                    sg_call = "Germline"
                
                funsig = get_funsig(af_dict, index)
                
                ## add new columns to row
                row['allAF'] = af_dict[index][0]
                row['Somatic/Germline Call'] = sg_call
                row['Functional significance'] = funsig

                ## write row to out_file
                writer.writerow(row)

main()
