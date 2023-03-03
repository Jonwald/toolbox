## parse TSO500 reports
import csv
import math
import pandas as pd
import sys
import os
from argparse import ArgumentParser

## samplesheet starts on 24

os.getcwd()

output_headers= ['Sample',
'Assay type',
'Variant type',
'Gene',
'SNV: Chromosome',
'SNV: Genomic Position',
'SNV: Reference Call',
'SNV: Alternative Call',
'SNV: Allele Frequency',
'SNV: Depth',
'SNV: P-Dot Notation',
'SNV: C-Dot Notation',
'SNV: Consequence(s)',
'SNV: Affected Exon(s)',
'CNV: Fold Change',
'TMB: Total TMB',
'MSI: Percent Unstable MSI Sites',
"FUSION: Gene Pair",
"FUSION: Breakpoint 1",
"FUSION: Breakpoint 2",
"FUSION: Fusion Supporting Reads",
"FUSION: Gene 1 Reference Reads",
"FUSION: Gene 2 Reference Reads",
"SPLICE: Affected Exon",
"SPLICE: Breakpoint 1",
"SPLICE: Breakpoint 2",
"SPLICE: Splice Supporting Reads",
"SPLICE: Reference Reads Transcript",
'DNA SNV/TMB QC: Median Exon Coverage (Count)',
'DNA SNV/TMB QC: Median Insert Size (bp)',
'DNA MSI QC: Usable MSI Sites (Count)',
'DNA SNV/TMB QC: % Exon 50X',
'DNA CNV QC: Coverage MAD Count',
'DNA CNV QC: Median Bin Count CNV Target (Count)',
"RNA QC: Median CV Gene 500X",
"RNA QC: Total On Target Reads",
"RNA QC: Median Insert Size"
]

def parse_command(command):
    """pasrses inputs"""
    parser = ArgumentParser(description='reformat TSO500 output for Bayer')
    parser.add_argument('-i', '--input', required=True, help='input folder path, containing "combinedvariantoutput and metricsoutput files')
    parser.add_argument('-o', '--output', required=True, help='output file prefix')

    try:
        args = parser.parse_args(command)
    except:
        parser.error('Invalid options provided')
    return (args)


def get_qc_table(csvfile, header, delim, type, assay, variant):
    with open(csvfile) as sample_sheet:
        header_index = False
        end_index = False
        readCSV = csv.reader(sample_sheet.readlines(), delimiter=delim)
        for line_num, content in enumerate(readCSV):
            if content and content[0] == header:
                header_index = line_num
            if header_index and not content[0]:
                end_index = line_num
                break
        header_index = header_index +1
        df = pd.read_csv(csvfile,
            skiprows= header_index,
            nrows=end_index - header_index -1,
            sep='\t')
    return df


#get index for given rows header form file
def get_table(csvfile, header, delim, type, assay, variant, sample):
    with open(csvfile) as sample_sheet:
        header_index = False
        end_index = False
        readCSV = csv.reader(sample_sheet.readlines(), delimiter=delim)
        for line_num, content in enumerate(readCSV):
            if content and content[0] == header:
                header_index = line_num
            if header_index and line_num > header_index and not content[1]: #and not content[0]:
                end_index = line_num
                break
        header_index = header_index +1
        #print(header_index,end_index)
        if type=="list":
            setnames=['Variable', 'value']
            df = pd.read_csv(csvfile,
            skiprows= header_index,
            nrows=end_index - header_index,
            sep='\t',
            names=setnames,
            index_col=False)
            df = df[df.columns.drop(list(df.filter(regex='Unnamed')))]
            df = pd.DataFrame(df.values.T[1:], columns=df.values.T[0])
        else:
            df = pd.read_csv(csvfile,
            skiprows= header_index,
            nrows=end_index - header_index -1,
            sep='\t')
        df = df.add_prefix(variant + ": ")
    df['Assay type'] = assay
    df['Variant type'] = variant
    df['Sample'] = sample
    return df

def convert_qc(table):
    dna_qc_call=[]
    table.columns = [c.replace(' ', '_') for c in table.columns]
    dna_cutoffs= table.LSL_Guideline.tolist()
    for index, row in table.iterrows():
        #print(index)
        tmp=[]
        #print(row[3:])
        for s in row[3:]:
            #print(s)
            if math.isnan(s):
                tmp.append("NA")
            elif s > dna_cutoffs[index]:
                tmp.append("PASS")
            else:
                tmp.append("FAIL")
        dna_qc_call.append(tmp)
    return dna_qc_call

def main():

    input_command = parse_command(sys.argv[1:])
    # get in data
    import glob
    sample_dat = glob.glob(os.path.join(input_command.input, "*_CombinedVariantOutput.tsv"))
    met_out = glob.glob(os.path.join(input_command.input,'*_MetricsOutput.tsv'))[0]
    #sample_dat = glob.glob("../misc/*_CombinedVariantOutput.tsv")
    #met_out = glob.glob("../misc/*_MetricsOutput.tsv")[0]
    sample_results=[]
    # loop over
    # SNV TMB MSI CNV
    # concatenate all results
    for sample in sample_dat:
        #sample="../misc/S2364-00019_CombinedVariantOutput.tsv"
        sname=os.path.basename(sample).replace("_CombinedVariantOutput.tsv", "")
        print(sname)
        dna_sv = get_table(sample, '[Small Variants]', '\t', "table", "DNA", "SNV", sname)
        dna_sv = dna_sv.rename(columns={ "SNV: Gene": "Gene" })
        dna_tmb = get_table(sample, '[TMB]', '\t', "list", "DNA", "TMB", sname)
        dna_msi = get_table(sample, '[MSI]', '\t',"list", "DNA", "MSI", sname)
        dna_cnv = get_table(sample, '[Gene Amplifications]', '\t', "table", "DNA", "CNV", sname)
        dna_cnv = dna_cnv.rename(columns={ "CNV: Gene": "Gene" })
        rna_fus = get_table(sample, '[Fusions]', '\t', "table", "RNA", "FUSION", sname)
        rna_spl = get_table(sample, '[Splice Variants]', '\t', "table", "RNA", "SPLICE", sname)
        ## add NA row to rna data if empty
        try:
            print(rna_fus.iloc[0])
        except:
            print("no fusions")
            rna_fus = rna_fus.append({"FUSION: Gene Pair": "NA",
            "FUSION: Breakpoint 1": "NA",
            "FUSION: Breakpoint 2": "NA",
            "FUSION: Fusion Supporting Reads": "NA",
            "FUSION: Gene 1 Reference Reads": "NA",
            "FUSION: Gene 2 Reference Reads": "NA"},
            ignore_index=True)

        try:
            print(rna_spl.iloc[0])
        except:
            print("no splice")
            rna_spl = rna_spl.append({"SPLICE: Gene": "NA",
            "SPLICE: Affected Exon": "NA",
            "SPLICE: Breakpoint 1": "NA",
            "SPLICE: Breakpoint 2": "NA",
            "SPLICE: Variant": "NA",
            "SPLICE: Splice Supporting Reads": "NA",
            "SPLICE: Reference Reads Transcript": "NA"},
            ignore_index=True)

        try:
            print(dna_cnv.iloc[0])
        except:
            print("no CNV")
            dna_cnv = dna_cnv.append({"Gene": "NA",
            "CNV: Fold Change": "NA",
            "CNV: Unnamed: 2": "NA",
            "Assay type": "NA",
            "Variant type": "NA",
            "Sample": "NA"},
            ignore_index=True)

        s_result_list=[dna_sv, dna_tmb, dna_msi, dna_cnv, rna_fus, rna_spl]
        sample_results.append(pd.concat(s_result_list, ignore_index=True))

    tidy_sample_results=pd.concat(sample_results, ignore_index=True)

    ## DNA QC
    qc_headers = ['[DNA Library QC Metrics for Small Variant Calling and TMB]',
    '[DNA Library QC Metrics for MSI]',
    '[DNA Library QC Metrics for CNV]']

    qc_tables = []
    for header in qc_headers:
        print(header)
        qc_tables.append(get_qc_table(met_out, header, '\t', "table", "DNA", "QC"))

    dna_qc = pd.concat(qc_tables, ignore_index=True)
    rna_qc = get_qc_table(met_out, '[RNA Library QC Metrics]', '\t', "table", "RNA", "QC")

    #dna_cutoffs[4] = dna_qc.USL_Guideline.tolist()[4]

    dna_qc_call = convert_qc(dna_qc)
    rna_qc_call = convert_qc(rna_qc)

    dna_headers=["DNA SNV/TMB QC: Median Insert Size (bp)",
    "DNA SNV/TMB QC: Median Exon Coverage (Count)",
    "DNA SNV/TMB QC: % Exon 50X",
    "DNA MSI QC: Usable MSI Sites (Count)",
    "DNA CNV QC: Coverage MAD Count",
    "DNA CNV QC: Median Bin Count CNV Target (Count)"]

    rna_headers=["RNA QC: Median CV Gene 500X",
    "RNA QC: Total On Target Reads",
    "RNA QC: Median Insert Size"]

    dna_p_f = pd.DataFrame(dna_qc_call)
    dna_tidy_qc = pd.DataFrame(dna_p_f.values.T, columns=dna_headers)
    dna_tidy_qc['Assay type'] = "DNA"
    dna_tidy_qc['Variant type'] = "QC"
    dna_tidy_qc['Sample'] = dna_qc.columns[3:].tolist()

    rna_p_f = pd.DataFrame(rna_qc_call)
    rna_tidy_qc = pd.DataFrame(rna_p_f.values.T, columns=rna_headers)
    rna_tidy_qc['Assay type'] = "RNA"
    rna_tidy_qc['Variant type'] = "QC"
    rna_tidy_qc['Sample'] = rna_qc.columns[3:].tolist()

    ## rna_qc

    final_out = pd.concat([dna_tidy_qc, rna_tidy_qc, tidy_sample_results], ignore_index=True)
    for col in final_out.columns:
        print(col)
    final_out = final_out.fillna("NA")

    final_out = final_out[((final_out.Sample != 'NA'))]
    final_out = final_out[((final_out.Sample != 'TSOCOMPDNA3'))]
    final_out = final_out[((final_out.Sample != 'TSOCOMPRNA1'))]
    final_out = final_out[((final_out.Sample != 'HD803'))]

    final_out[output_headers].to_csv(input_command.output + ".tsv", sep="\t", index=False)

main()

