import os 
import vcf
import csv

"""
insert sample ids here
sample_ids = ["S1910-00016", "S1910-00027", "S1984-00004", "S1984-00012", "S1984-00021", "S1984-00024", "S1984-00026", "S1984-00029",
"S1984-00047", "S1984-00048", "S1984-00050", "S1984-00066", "S1984-00071", "S1984-00073", "S1984-00074", "S1984-00076",
"S1984-00087", "S1984-00089", "S1984-00097", "S1984-00098", "S1984-00099"]
"""

header = ["sample ID",
"Chromosome",
"Position",
"Reference",
"Allele",
"Allele type",  #(indel, snv, …)
"cDNA change",
"Protein change", # if any
"Variant allele frequency",
"Allele Read depth",
"Dbsnp id",  #if any
"COSMIC id",  # if any
"Clinvar annotation", # (ClinVar significance, benign, etc) if any and 
"EXAC freq"] # if any

with open("CYT_somatic_variants.csv", 'w') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for sample_id in sample_ids:
        in_file = sample_id + '.somatic.filtered.vcf'
        print("processing: " + sample_id)
        vcf_reader = vcf.Reader(open(in_file, 'r'))    
        for rec in vcf_reader:
            #print(rec.POS, rec.var_type)
            if 'vep' in rec.INFO.keys():
                cdna = rec.INFO['vep'][0].split('|')[10]
                protein = rec.INFO['vep'][0].split('|')[11]
                protein_id = rec.INFO['vep'][0].split('|')[3]
            else:
                cdna = "NA"
                protein = "NA"
                protein_id = "NA"
            if cdna =='':
                cdna= "NA"
            if  protein =='':
                protein = "NA"     
            for sample in rec.samples:
                AF = sample['AF']
                DP = sample['DP']
            if "RS" in rec.INFO.keys():
                DBSNP_tmp = rec.INFO['RS']
                DBSNP_ID = []
                for i in DBSNP_tmp:
                    DBSNP_ID.append('rs' + i)
            else:
                DBSNP_ID = "NA"
            if len(DBSNP_ID) > 1:
                DBSNP_ID = '|'.join(DBSNP_ID)
            else:
                DBSNP_ID = DBSNP_ID[0]
            if "LEGACY_ID" in rec.INFO.keys():
                COSMIC_ID = rec.INFO['LEGACY_ID']
            else:
                COSMIC_ID = "NA"
            if "CLNSIG" in rec.INFO.keys():
                CLINSIG = rec.INFO['CLNSIG']
            else:
                CLINSIG = "NA"
            if len(CLINSIG) > 1:
                CLINSIG = '|'.join(CLINSIG)
            else:
                CLINSIG = CLINSIG[0]
            if "AF" in rec.INFO.keys():
                EXAC_FREQ= rec.INFO['AF']
            else:
                EXAC_FREQ = "NA"
            if len(EXAC_FREQ) > 1:
                EXAC_FREQ = '|'.join(str(EXAC_FREQ))
            else:
                EXAC_FREQ = EXAC_FREQ[0]
            if len(rec.ALT) > 1:
                ALT_tmp = []
                for i in rec.ALT:
                    ALT_tmp.append(str(i))
                ALT = '|'.join(ALT_tmp)
            else:
                ALT = rec.ALT[0]
            line = [sample_id, rec.CHROM, rec.POS, rec.REF, ALT, rec.var_type, cdna, protein, AF, DP, DBSNP_ID, COSMIC_ID, CLINSIG, EXAC_FREQ]
            writer.writerow(line)


## germline files
with open("CYT_germline_variants.csv", 'w') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for sample_id in sample_ids:
        in_file = sample_id + '.germline.filtered.vcf'
        print("processing: " + sample_id)
        vcf_reader = vcf.Reader(open(in_file, 'r'))
        for rec in vcf_reader:
            #print(rec.POS, rec.var_type)
            if 'vep' in rec.INFO.keys():
                cdna = rec.INFO['vep'][0].split('|')[10]
                protein = rec.INFO['vep'][0].split('|')[11]
            else:
                cdna = "NA"
                protein = "NA"
            if cdna =='':
                cdna= "NA"
            if  protein =='':
                protein = "NA"     
            for s in rec.samples:
                DP = s['DP']
                AD = s['AD']
                AD_tmp = []
                for k in test:
                    AD_tmp.append(str(round(k/sum(s['AD']), 3)))
                AF = AD_tmp[1:]
                AF = '|'.join(AF)
            if "RS" in rec.INFO.keys():
                DBSNP_tmp = rec.INFO['RS']
                DBSNP_ID = []
                for i in DBSNP_tmp:
                    DBSNP_ID.append('rs' + i)
            else:
                DBSNP_ID = "NA"
            if len(DBSNP_ID) > 1:
                DBSNP_ID = '|'.join(DBSNP_ID)
            else:
                DBSNP_ID = DBSNP_ID[0]
            if "LEGACY_ID" in rec.INFO.keys():
                COSMIC_ID = rec.INFO['LEGACY_ID']
            else:
                COSMIC_ID = "NA"
            if "CLNSIG" in rec.INFO.keys():
                CLINSIG = rec.INFO['CLNSIG']
            else:
                CLINSIG = "NA"
            if len(CLINSIG) > 1:
                CLINSIG = '|'.join(CLINSIG)
            else:
                CLINSIG = CLINSIG[0]
            if "AN" in rec.INFO.keys():
                EXAC_FREQ= rec.INFO['AF']
            else:
                EXAC_FREQ = "NA"
            if len(EXAC_FREQ) > 1:
                EXAC_FREQ = '|'.join(str(EXAC_FREQ))
            else:
                EXAC_FREQ = EXAC_FREQ[0]
            if len(rec.ALT) > 1:
                ALT_tmp = []
                for i in rec.ALT:
                    ALT_tmp.append(str(i))
                ALT = '|'.join(ALT_tmp)
            else:
                ALT = rec.ALT[0]
            line = [sample_id, rec.CHROM, rec.POS, rec.REF, ALT, rec.var_type, cdna, protein, AF, DP, DBSNP_ID, COSMIC_ID, CLINSIG, EXAC_FREQ]
            writer.writerow(line)



test = [0, 5, 10]
test_out = []
for i in test:
    test_out.append(round(i/sum(test), 3))



import os 
import vcf
import csv
os.getcwd()


sample_ids = ["S1910-00016", "S1910-00027", "S1984-00004", "S1984-00012", "S1984-00021", "S1984-00024", "S1984-00026", "S1984-00029",
"S1984-00047", "S1984-00048", "S1984-00050", "S1984-00066", "S1984-00071", "S1984-00073", "S1984-00074", "S1984-00076",
"S1984-00087", "S1984-00089", "S1984-00097", "S1984-00098", "S1984-00099"]

for sample_id in sample_ids:
    sample_id = sample_ids[0]
    print("processing sample: " + sample_id)
    germline_file = sample_id + ".germline.filtered.vcf"
    vcf_reader = vcf.Reader(open(germline_file, 'r'))  
    #vcf_writer = vcf.Writer(open(sample_id + '_1p_filtered.vcf', 'w'), vcf_reader)
    for rec in vcf_reader:
        #vcf_headers = rec.INFO.keys()
        if 'vep' in rec.INFO.keys():
            protein = rec.INFO['vep']
            protein[0].split('|')[3]

        vcf_headers
        for s in rec.samples:
            
            DP = s['DP']
            AD_tmp = []
            for k in test:
                AD_tmp.append(round(k/sum(s['AD']), 3))
            AF = AD_tmp[1:]
            if any([x > 0.01 for x in AF]):
                vcf_writer.write_record(record)



