#!/usr/bin/env python

"""
This is the testing unit for FishTaco
"""
# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import unittest
import subprocess
import os
import glob

# when the test module will be ready:
import fishtaco
from fishtaco.compute_contribution_to_DA import main

# for testing:
#import sys
#sys.path.append('/net/gs/vol1/home/ohadm/METAFIT/PyCode/FiShTaCo/fishtaco')
#from compute_contribution_to_DA import main


class FishTacoTestCase(unittest.TestCase):
    """Tests for `compute_contribution_to_DA.py`."""

    # Get the path to the location of the package:
    path_to_data = os.path.dirname(fishtaco.__file__)

    # for testing:
    # path_to_data = '/net/gs/vol1/home/ohadm/METAFIT/PyCode/FiShTaCo/fishtaco'

    def test_is_output_correct_for_fishtaco_no_inference(self):
        """Does FishTaco with no inference produce the correct output for the example case?"""
        print("Path to examples:" + FishTacoTestCase.path_to_data)

        # run fishtaco from shell to make sure it works properly
        subprocess.call('run_fishtaco.py -op fishtaco_out_no_inf -mult_hyp Bonf' +
                        ' -ta ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data + '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data + '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter K00001 ' +
                        ' -map_function_level none -functional_profile_already_corrected_with_musicc', shell=True)

        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_no_inf_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')

        print("Testing output...")
        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_no_inf_*.tab"):
            os.remove(fl)

        self.assertTrue('Program completed successfully' in output)

    def test_is_output_correct_for_fishtaco_prior_based_inference(self):
        """Does FishTaco with prior-based inference produce the correct output for the example case?"""

        # run fishtaco from shell to make sure it works properly
        subprocess.call('run_fishtaco.py -op fishtaco_out_prior_based_inf -mult_hyp Bonf -inf' +
                        ' -ta ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data + '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data + '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter K00001 ' +
                        ' -map_function_level none -functional_profile_already_corrected_with_musicc', shell=True)

        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_prior_based_inf_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')

        print("Testing output...")
        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_prior_based_inf_*.tab"):
            os.remove(fl)

        self.assertTrue('Program completed successfully' in output)

    def test_is_output_correct_for_fishtaco_de_novo_inference(self):
        """Does FishTaco with de novo inference produce the correct output for the example case?"""

        # run fishtaco from shell to make sure it works properly
        subprocess.call('run_fishtaco.py -op fishtaco_out_de_novo_inf -mult_hyp Bonf -inf' +
                        ' -ta ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data + '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data + '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter K00001 ' +
                        ' -map_function_level none -functional_profile_already_corrected_with_musicc', shell=True)

        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_de_novo_inf_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')

        print("Testing output...")
        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_de_novo_inf_*.tab"):
            os.remove(fl)

        self.assertTrue('Program completed successfully' in output)

    def test_is_output_correct_for_fishtaco_predict_functional_profile(self):
        """Does FishTaco with predicting the functional profiles produce the correct output for the example case?"""

        # run fishtaco with the option to predict function (by not supplying a functional profile to the script)
        subprocess.call('run_fishtaco.py -op fishtaco_out_predict_function -mult_hyp Bonf' +
                        ' -ta ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data + '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter K00001 ' +
                        ' -map_function_level none', shell=True)

        # run fishtaco with inference on the functional profile we have just predicted
        subprocess.call('run_fishtaco.py -op fishtaco_out_infer_from_predicted_function -mult_hyp Bonf -inf' +
                        ' -ta ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu fishtaco_out_predict_function_STAT_predicted_function_abundance_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data + '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter K00001 ' +
                        ' -map_function_level none -functional_profile_already_corrected_with_musicc', shell=True)

        # assert that the log files show that fishtaco completed successfully
        with open("fishtaco_out_predict_function_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')
        self.assertTrue('Program completed successfully' in output)

        with open("fishtaco_out_infer_from_predicted_function_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')
        self.assertTrue('Program completed successfully' in output)

        print("Testing output...")
        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_predict_function_*.tab"):
            os.remove(fl)
        for fl in glob.glob("fishtaco_out_infer_from_predicted_function_*.tab"):
            os.remove(fl)

    def test_is_output_correct_for_fishtaco_filtering_by_function_list(self):
        """Does FishTaco with filtering by function list produce the correct output for the example case?"""
        print("Path to examples:" + FishTacoTestCase.path_to_data)

         # run fishtaco without filtering
        subprocess.call('run_fishtaco.py -op fishtaco_out_no_filtering -mult_hyp Bonf' +
                        ' -ta ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data + '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001_K00054.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data + '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_KO_only_K00001_K00054.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5' +
                        ' -map_function_level none -functional_profile_already_corrected_with_musicc', shell=True)

        # run fishtaco with filtering a list of functions
        subprocess.call('run_fishtaco.py -op fishtaco_out_filtering_by_function_list -mult_hyp Bonf' +
                        ' -ta ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data + '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001_K00054.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data + '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data + '/examples/METAPHLAN_taxa_vs_KO_only_K00001_K00054.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -multi_function_filter_list K00001,K00007,K00020' +
                        ' -map_function_level none -functional_profile_already_corrected_with_musicc', shell=True)

        print("Testing output...")
        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_no_filtering_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')
        self.assertTrue('Program completed successfully' in output)

        with open("fishtaco_out_filtering_by_function_list_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')
        self.assertTrue('Program completed successfully' in output)

        # count number of lines in the filtered file and make sure it is equal to
        with open('fishtaco_out_filtering_by_function_list_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab') as my_file:
            num_lines = len(my_file.readlines())
        self.assertTrue(num_lines == 8)

        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_no_filtering_*.tab"):
            os.remove(fl)
        for fl in glob.glob("fishtaco_out_filtering_by_function_list_*.tab"):
            os.remove(fl)

################################################

if __name__ == '__main__':
    unittest.main()


