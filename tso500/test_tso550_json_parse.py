import unittest
import csv
import io
import sys
import os
import json
import gzip
from argparse import ArgumentParser
import collections
from ts500_json_parse import get_file_pairs
from ts500_json_parse import read_json
from ts500_json_parse import parse_tmb
from ts500_json_parse import parse_command
from ts500_json_parse import get_sv_start
from ts500_json_parse import copy_first_x_lines
from ts500_json_parse import is_somatic


class TestParseCommand(unittest.TestCase):
    def test_parse_command(self):
        # Test the parse_command function with valid arguments
        command = ['--input', '/path/to/input', '--output', '/path/to/output']
        args = parse_command(command)
        self.assertEqual(args.input, '/path/to/input')
        self.assertEqual(args.output, '/path/to/output')

        # Test the parse_command function with missing arguments
        command = ['--input', '/path/to/input']
        with self.assertRaises(SystemExit):
            parse_command(command)


class TestParseTMB(unittest.TestCase):
    def setUp(self):
        # Create a mock TMB file for testing
        self.mock_tmb_file = 'mock_tmb_file.tsv'
        with open(self.mock_tmb_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['Chromosome', 'Position', 'RefCall', 'AltCall', 'MaxCosmicCount'])
            writer.writerow(['chr1', '12345', 'A', 'T', '10'])
            writer.writerow(['chr2', '67890', 'C', 'G', '20'])

    def tearDown(self):
        # Remove the mock TMB file after testing
        os.remove(self.mock_tmb_file)

    def test_parse_tmb(self):
        # Test the parse_tmb function
        expected_tmb_dict = {'chr1_12345_A_T': 10,
                             'chr2_67890_C_G': 20}
        self.assertEqual(parse_tmb(self.mock_tmb_file), expected_tmb_dict)


class TestReadJson(unittest.TestCase):
    def setUp(self):
        # Create a mock JSON file for testing
        self.mock_json_file = 'mock_json_file.json.gz'
        data = {
            # pass variant with gnomad and clinvar annpotations
            'positions': [
                {'filters': ['PASS'],
                'chromosome': 'chr1',
                'position': 12345,
                'refAllele': 'A',
                'altAlleles': ['T'],
                'variants': [
                            {'vid': '1:12345:T',
                            'gnomad': {'allAc': 10, 'allAn': 20},
                            'gnomadExome': {'allAc': 5, 'allAn': 10},
                            'clinvar': [{'alleleOrigins': ['inherited', 'maternal']}]}
                            ],
                
                },
            # pass variant with no gnomad or clinvar annotations
                {'filters': ['PASS'],
                'chromosome': 'chr2',
                'position': 67890,
                'refAllele': 'C',
                'altAlleles': ['G'],
                'variants': [{'vid': '2:67890:G'}]
                },
            # fail variant
                {'filters': ['LowSupport'],
                'chromosome': 'chr3',
                'position': 12345,
                'refAllele': 'T',
                'altAlleles': ['A'],
                'variants': [{'vid': '3:12345:A'}]}
            ]
        }

        with gzip.open(self.mock_json_file, 'wt') as f:
            json.dump(data, f)
        
    def tearDown(self):
        # Remove the mock JSON file after testing
        os.remove(self.mock_json_file)

    def test_read_json(self):
        # Test the read_json function
        tmb_dict = {'chr1_12345_A_T': 10, 'chr2_67890_C_G': 0}
        expected_af_dict = {'chr1_12345_A_T': [0.5, 10, ['inherited', 'maternal']],
                            'chr2_67890_C_G': [0, 0, ['NA']]}
        self.assertEqual(read_json(self.mock_json_file, tmb_dict), expected_af_dict)


class TestGetSVStart(unittest.TestCase):
    def setUp(self):
        # Create a mock file for testing
        self.mock_file = 'mock_file.txt'
        with open(self.mock_file, 'w') as f:
            f.write('Some text\n')
            f.write('More text\n')
            f.write('[Small Variants]\n')
            f.write('Yet more text\n')

    def tearDown(self):
        # Remove the mock file after testing
        os.remove(self.mock_file)

    def test_get_sv_start(self):
        # Test the get_sv_start function
        self.assertEqual(get_sv_start(self.mock_file), 3)

class TestIsSomatic(unittest.TestCase):
    def test_is_somatic(self):
        # consequence filter
        cons_filter = ['missense_variant', 'frameshift_variant', 'splice_acceptor_variant',
                     'splice_donor_variant', 'inframe_deletion','inframe_insertion',
                     'stop_gained', 'stop_lost', 'start_lost', '3_prime_UTR_variant',
                     '5_prime_UTR_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                     'transcript_ablation', 'transcript_amplification', 'feature_elongation',
                      'feature_truncation', 'protein_altering_variant']

        # clinvar germline annotations
        clin_ann = ['inherited', 'paternal', 'maternal', 'biparental', 'uniparental']

        # Test 'True' case 1 
        # match to funn_ann, allAF < 0.01, cosmic count > 20
        index = 'chr1_12345_A_T'
        af_dict = {'chr1_12345_A_T': [0.005, 25, ['somatic']]}
        fun_ann = ['frameshift_variant']
        vf = 0.8
        self.assertTrue(is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf))

        # Test 'True' case 2
        # match to funn_ann, allAF < 0.01, cosmic count < 20, but VF < 0.9 and no clin_ann match
        fun_ann = ['frameshift_variant']
        vf = 0.8
        af_dict = {'chr1_12345_A_T': [0.005, 25, ['unknown']]}
        self.assertTrue(is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf))
        
        # Test 'False' case 1 
        # allAF not match to fun_ann
        fun_ann = ['intron_variant']
        vf = 0.8
        af_dict = {'chr1_12345_A_T': [0.02, 25, ['unknown']]}
        self.assertFalse(is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf))

        # Test 'False' case 2 
        # allAF > 0.01
        fun_ann = ['frameshift_variant']
        vf = 0.8
        af_dict = {'chr1_12345_A_T': [0.02, 25, ['unknown']]}
        self.assertFalse(is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf))

        # Test 'False' case 3 
        # cosmic count < 20, vf > 0.9
        fun_ann = ['frameshift_variant']
        vf = 0.95
        af_dict = {'chr1_12345_A_T': [0.005, 15, ['unknown']]}
        self.assertFalse(is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf))

        # Test 'False' case 4
        # cosmic count < 20, vf < 0.9 but clin_ann match
        fun_ann = ['frameshift_variant']
        vf = 0.8
        af_dict = {'chr1_12345_A_T': [0.005, 15, ['inherited']]}
        self.assertFalse(is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf))

        # Test 'False' case 5
        # multple clinvar annotations, one matches to clin_ann
        af_dict = {'chr1_12345_A_T': [0.005, 15, ['somatic', 'somatic', 'inherited']]}
        vf = 0.8
        self.assertFalse(is_somatic(fun_ann, af_dict, index, cons_filter, clin_ann, vf))        

class TestCopyFirstXLines(unittest.TestCase):
    def setUp(self):
        # Create a mock input file for testing
        self.mock_input_file = 'mock_input_file.tsv'
        with open(self.mock_input_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            for i in range(10):
                writer.writerow([i, 'row' + str(i)])

        # Create a mock output file for testing
        self.mock_output_file = 'mock_output_file.tsv'

    def tearDown(self):
        # Remove the mock files after testing
        os.remove(self.mock_input_file)
        os.remove(self.mock_output_file)

    def test_copy_first_x_lines(self):
        # Test the copy_first_x_lines function
        copy_first_x_lines(self.mock_input_file, self.mock_output_file, 5)
        with open(self.mock_output_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)
            self.assertEqual(len(rows), 5)
            for i, row in enumerate(rows):
                self.assertEqual(row, [str(i), 'row' + str(i)])

class TestGetFilePairs(unittest.TestCase):
    def setUp(self):
        # Create a mock directory for testing
        self.mock_directory = 'mock_directory'
        os.mkdir(self.mock_directory)

        # Create some mock files in the directory
        self.mock_files = ['sample1_MergedVariants_Annotated.json.gz',
                           'sample1_CombinedVariantOutput.tsv',
                           'sample1_TMB_Trace.tsv',
                           'sample2_MergedVariants_Annotated.json.gz',
                           'sample2_CombinedVariantOutput.tsv',
                           'sample2_TMB_Trace.tsv',
                           'sample3_MergedVariants_Annotated.json.gz',
                           'sample3_CombinedVariantOutput.tsv',
                           'sample3_TMB_Trace.tsv']
        for mock_file in self.mock_files:
            with open(os.path.join(self.mock_directory, mock_file), 'w') as f:
                f.write('')

    def tearDown(self):
        # Remove the mock directory and its contents after testing
        for mock_file in self.mock_files:
            os.remove(os.path.join(self.mock_directory, mock_file))
        os.rmdir(self.mock_directory)

    def test_get_file_pairs(self):
        # Test the get_file_pairs function
        expected_result = collections.defaultdict(list)
        expected_result['sample1'] = set(['sample1_CombinedVariantOutput.tsv',
                                      'sample1_MergedVariants_Annotated.json.gz',
                                      'sample1_TMB_Trace.tsv'])
        expected_result['sample2'] = set(['sample2_CombinedVariantOutput.tsv',
                                      'sample2_MergedVariants_Annotated.json.gz',
                                      'sample2_TMB_Trace.tsv'])
        expected_result['sample3'] = set(['sample3_CombinedVariantOutput.tsv',
                                      'sample3_MergedVariants_Annotated.json.gz',
                                      'sample3_TMB_Trace.tsv'])

        self.assertDictEqual(get_file_pairs(self.mock_directory), expected_result)

if __name__ == '__main__':
    unittest.main(verbosity=2)

