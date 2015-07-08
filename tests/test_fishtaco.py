#!/usr/bin/env python

"""
This is the testing unit for FishTaco
"""
# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import unittest
import subprocess
import os

# when the test module will be ready:
import fishtaco
from fishtaco.compute_contribution_to_DA import main

# for testing:
#import sys
#sys.path.append('/net/gs/vol1/home/ohadm/METAFIT/PyCode/FiShTaCo/fishtaco')
#from compute_contribution_to_DA import main


class FishTacoTestCase(unittest.TestCase):
    """Tests for `compute_contribution_to_DA.py`."""

    # when the test module will be ready:
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
        os.remove("fishtaco_out_no_inf_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab")

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
        os.remove("fishtaco_out_prior_based_inf_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab")

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
        os.remove("fishtaco_out_de_novo_inf_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_single_taxa.tab")

        self.assertTrue('Program completed successfully' in output)

################################################

if __name__ == '__main__':
    unittest.main()


