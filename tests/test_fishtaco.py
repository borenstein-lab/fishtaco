#!/usr/bin/env python

"""
This is the testing unit for FishTaco
"""
# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import unittest
import subprocess
import os
import sys
import glob
import pandas as pd
import numpy as np
import fishtaco


class FishTacoTestCase(unittest.TestCase):
    """Tests for `compute_contribution_to_DA.py`.
    """

    # Get the path to the location of the package:
    path_to_data = os.path.dirname(fishtaco.__file__)

    def test_compute_differential_abundance(self):
        """Does compute_differential_abundance.py produce the correct output
        for the example case?"""
        print("==============================================================")
        print("Testing compute_differential_abundance.py")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run compute_differential_abundance from shell to make sure it
        # works properly
        subprocess.call('python ' +
                        FishTacoTestCase.path_to_data +
                        '/compute_differential_abundance.py ' +
                        FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -c ' +
                        FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab -ch ' +
                        '-o test_compute_differential_abundance.tab -v',
                        shell=True)

        example_output = \
            pd.read_csv(FishTacoTestCase.path_to_data +
                        '/examples/output/fishtaco_out_no_inf_STAT_DA_taxa_'
                        'SCORE_wilcoxon_ASSESSMENT_single_taxa.tab', sep="\t",
                        comment="#")

        test_output = \
            pd.read_csv('test_compute_differential_abundance.tab', sep="\t")

        example_stat_values = example_output[["StatValue"]].values.flatten()
        test_stat_values = \
            test_output.sort_values(by="Function")[["StatValue"]].values.flatten()

        os.remove('test_compute_differential_abundance.tab')

        self.assertTrue(np.allclose(example_stat_values, test_stat_values))

    def test_compute_pathway_abundance(self):
        """Does compute_pathway_abundance.py produce the correct output
        for the example case?"""
        print("==============================================================")
        print("Testing compute_pathway_abundance.py")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run compute_pathway_abundance from shell to make sure it
        # works properly
        subprocess.call('python ' +
                        FishTacoTestCase.path_to_data +
                        '/compute_pathway_abundance.py ' +
                        ' -ko ' + FishTacoTestCase.path_to_data +
                        '/examples/KO_vs_SAMPLE_MUSiCC.tab' +
                        ' -ko2path ' +
                        FishTacoTestCase.path_to_data +
                        '/data/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15.tab ' +
                        '-o test_compute_pathway_abundance.tab ' +
                        '-oc test_compute_pathway_abundance_counts.tab -v',
                        shell=True)

        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                        '/examples/PATHWAY_vs_SAMPLE_MUSiCC.tab', sep="\t")

        test_output = pd.read_csv('test_compute_pathway_abundance.tab',
                                  sep="\t")

        example_values = example_output.values[:, 1:5].flatten().astype(float)
        test_values = test_output.values[:, 1:5].flatten().astype(float)

        os.remove('test_compute_pathway_abundance.tab')
        os.remove('test_compute_pathway_abundance_counts.tab')

        self.assertTrue(np.allclose(example_values, test_values))

    def test_learn_non_neg_elastic_net_with_prior(self):
        """Does learn_non_neg_elastic_net_with_prior.py produce the correct
        output for the example case when we use a prior?"""
        print("==============================================================")
        print("Testing learn_non_neg_elastic_net_with_prior.py (prior-based)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run FishTaco with inference from shell and examine that the
        # inference was done correctly. First use prior in inference:
        subprocess.call('run_fishtaco.py -op test_prior_based_inf '
                        '-mult_hyp Bonf -inf' +
                        ' -ta ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data +
                        '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter '
                        'K00001 ' +
                        ' -map_function_level none '
                        '-functional_profile_already_corrected_with_musicc',
                        shell=True)

        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_prior_'
                                     'based_inf_STAT_taxa_learned_copy_num_'
                                     'SCORE_wilcoxon_ASSESSMENT_single_'
                                     'taxa.tab', sep="\t")

        test_output = pd.read_csv('test_prior_based_inf_STAT_taxa_learned_'
                                  'copy_num_SCORE_wilcoxon_ASSESSMENT_single_'
                                  'taxa.tab', sep="\t")

        example_values = example_output.values[:, 1].flatten().astype(float)
        test_values = test_output.values[:, 1].flatten().astype(float)

        for fl in glob.glob("test_prior_based_inf*.tab"):
            os.remove(fl)

        self.assertTrue(np.allclose(example_values, test_values))

    def test_learn_non_neg_elastic_net_de_novo(self):
        """Does learn_non_neg_elastic_net_with_prior.py produce the correct
        output for the example case when we DO NOT use a prior (de novo)?"""
        print("==============================================================")
        print("Testing learn_non_neg_elastic_net_with_prior.py (de-novo)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run FishTaco with inference from shell and examine that the
        # inference was done correctly. No prior in inference:
        subprocess.call('run_fishtaco.py -op test_de_novo_based_inf '
                        '-mult_hyp Bonf -inf' +
                        ' -ta ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data +
                        '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter '
                        'K00001 ' +
                        ' -map_function_level none '
                        '-functional_profile_already_corrected_with_musicc',
                        shell=True)

        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_de_novo_'
                                     'inf_STAT_taxa_learned_copy_num_'
                                     'SCORE_wilcoxon_ASSESSMENT_single_'
                                     'taxa.tab', sep="\t")

        test_output = pd.read_csv('test_de_novo_based_inf_STAT_taxa_learned_'
                                  'copy_num_SCORE_wilcoxon_ASSESSMENT_single_'
                                  'taxa.tab', sep="\t")

        example_values = example_output.values[:, 1].flatten().astype(float)
        test_values = test_output.values[:, 1].flatten().astype(float)

        for fl in glob.glob("test_de_novo_based_inf*.tab"):
            os.remove(fl)

        self.assertTrue(np.allclose(example_values, test_values))

    def test_is_output_correct_for_fishtaco_no_inference(self):
        """Does FishTaco with no inference produce the correct output for
        the example case?
        """
        print("==============================================================")
        print("Testing compute_contribution_to_DA.py (no inference)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run fishtaco from shell to make sure it works properly
        subprocess.call('run_fishtaco.py -op fishtaco_out_no_inf -mult_hyp '
                        'Bonf' +
                        ' -ta ' +
                        FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data +
                        '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter '
                        'K00001 ' +
                        ' -map_function_level none '
                        '-functional_profile_already_corrected_with_musicc',
                        shell=True)

        print("Testing output...")

        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_no_inf_STAT_run_log_SCORE_wilcoxon"
                  "_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')

        self.assertTrue('Program completed successfully' in output)

        # assert that the calculations of the shift contribution of each taxon
        # are the same as the example
        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_no_inf_'
                                     'STAT_taxa_contributions_SCORE_wilcoxon_'
                                     'ASSESSMENT_single_taxa.tab', sep="\t")

        test_output = pd.read_csv('fishtaco_out_no_inf_'
                                  'STAT_taxa_contributions_SCORE_wilcoxon_'
                                  'ASSESSMENT_single_taxa.tab', sep="\t")

        example_values = example_output.values[:, 1].flatten().astype(float)
        test_values = test_output.values[:, 1].flatten().astype(float)

        self.assertTrue(np.allclose(example_values, test_values, atol=2.0))

        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_no_inf_*.tab"):
            os.remove(fl)

    def test_is_output_correct_for_fishtaco_prior_based_inference(self):
        """Does FishTaco with prior-based inference produce the correct
        output for the example case?
        """
        print("==============================================================")
        print("Testing compute_contribution_to_DA.py (prior-based inference)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run fishtaco from shell to make sure it works properly
        subprocess.call('run_fishtaco.py -op fishtaco_out_prior_based_inf '
                        '-mult_hyp Bonf -inf' +
                        ' -ta ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data +
                        '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter '
                        'K00001 ' +
                        ' -map_function_level none '
                        '-functional_profile_already_corrected_with_musicc',
                        shell=True)

        print("Testing output...")

        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_prior_based_inf_STAT_run_log_SCORE_wilcoxon"
                  "_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')

        self.assertTrue('Program completed successfully' in output)

        # assert that the calculations of the shift contribution of each taxon
        # are the same as the example
        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_prior_'
                                     'based_inf_'
                                     'STAT_taxa_contributions_SCORE_wilcoxon_'
                                     'ASSESSMENT_single_taxa.tab', sep="\t")

        test_output = pd.read_csv('fishtaco_out_prior_based_inf_'
                                  'STAT_taxa_contributions_SCORE_wilcoxon_'
                                  'ASSESSMENT_single_taxa.tab', sep="\t")

        example_values = example_output.values[:, 1].flatten().astype(float)
        test_values = test_output.values[:, 1].flatten().astype(float)

        self.assertTrue(np.allclose(example_values, test_values, atol=2.0))

        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_prior_based_inf_*.tab"):
            os.remove(fl)

    def test_is_output_correct_for_fishtaco_de_novo_inference(self):
        """Does FishTaco with de novo inference produce the correct output
        for the example case?
        """
        print("==============================================================")
        print("Testing compute_contribution_to_DA.py (de novo inference)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run fishtaco from shell to make sure it works properly
        subprocess.call('run_fishtaco.py -op fishtaco_out_de_novo_inf '
                        '-mult_hyp Bonf -inf' +
                        ' -ta ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data +
                        '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter '
                        'K00001 ' +
                        ' -map_function_level none  '
                        '-functional_profile_already_corrected_with_musicc',
                        shell=True)

        print("Testing output...")

        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_de_novo_inf_STAT_run_log_SCORE_wilcoxon"
                  "_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')

        self.assertTrue('Program completed successfully' in output)

        # assert that the calculations of the shift contribution of each taxon
        # are the same as the example
        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_de_novo_'
                                     'inf_STAT_taxa_contributions_SCORE_'
                                     'wilcoxon_'
                                     'ASSESSMENT_single_taxa.tab', sep="\t")

        test_output = pd.read_csv('fishtaco_out_de_novo_inf_'
                                  'STAT_taxa_contributions_SCORE_wilcoxon_'
                                  'ASSESSMENT_single_taxa.tab', sep="\t")

        example_values = example_output.values[:, 1].flatten().astype(float)
        test_values = test_output.values[:, 1].flatten().astype(float)

        self.assertTrue(np.allclose(example_values, test_values, atol=2.0))

        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_de_novo_inf_*.tab"):
            os.remove(fl)

    def test_is_output_correct_for_fishtaco_predict_functional_profile(self):
        """Does FishTaco with predicting the functional profiles produce the
        correct output for the example case?
        """
        print("==============================================================")
        print("Testing compute_contribution_to_DA.py (predict function)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run fishtaco with the option to predict function (by not supplying a
        # functional profile to the script)
        subprocess.call('run_fishtaco.py -op fishtaco_out_predict_function '
                        '-mult_hyp Bonf' +
                        ' -ta ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 -single_function_filter '
                        'K00001 ' +
                        ' -map_function_level none', shell=True)

        print("Testing output...")

        # assert that the log files show that fishtaco completed successfully
        with open("fishtaco_out_predict_function_STAT_run_log_SCORE_wilcoxon"
                  "_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')
        self.assertTrue('Program completed successfully' in output)

        # assert that the calculations of the shift contribution of each taxon
        # are the same as the example
        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_predict_'
                                     'function_STAT_predicted_function_'
                                     'abundance_SCORE_wilcoxon_'
                                     'ASSESSMENT_single_taxa.tab', sep="\t")

        test_output = pd.read_csv('fishtaco_out_predict_function_'
                                  'STAT_predicted_function_abundance_SCORE_'
                                  'wilcoxon_ASSESSMENT_single_taxa.tab',
                                  sep="\t")

        example_values = example_output.values[0, 1:].flatten().astype(float)
        test_values = test_output.values[0, 1:].flatten().astype(float)

        self.assertTrue(np.allclose(example_values, test_values, atol=2.0))

        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_predict_function_*.tab"):
            os.remove(fl)

    def test_is_output_correct_for_fishtaco_filtering_by_function_list(self):
        """Does FishTaco with filtering by function list produce the correct
        output for the example case?
        """
        print("==============================================================")
        print("Testing compute_contribution_to_DA.py (filter by list)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run fishtaco with filtering a list of functions
        subprocess.call('run_fishtaco.py -op '
                        'fishtaco_out_filtering_by_function_list -mult_hyp Bonf' +
                        ' -ta ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data +
                        '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001_K00054.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_KO_only_K00001_K00054.tab' +
                        ' -assessment single_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 100 -log ' +
                        ' -number_of_permutations 5 '
                        '-multi_function_filter_list K00001,K00007,K00020' +
                        ' -map_function_level none '
                        '-functional_profile_already_corrected_with_musicc',
                        shell=True)

        print("Testing output...")
        # assert that the log file shows that fishtaco completed successfully
        with open("fishtaco_out_filtering_by_function_list_STAT_run_log_SCORE"
                  "_wilcoxon_ASSESSMENT_single_taxa.tab", "r") as my_file:
            output = my_file.read().replace('\n', '')
        self.assertTrue('Program completed successfully' in output)

        # assert that the calculations of the shift contribution of each taxon
        # are the same as the example
        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_filtering_'
                                     'by_function_list_STAT_taxa_contributions_'
                                     'SCORE_wilcoxon_'
                                     'ASSESSMENT_single_taxa.tab', sep="\t")

        test_output = pd.read_csv('fishtaco_out_filtering_by_function_list_'
                                  'STAT_taxa_contributions_SCORE_'
                                  'wilcoxon_ASSESSMENT_single_taxa.tab',
                                  sep="\t")

        example_values = example_output.values[:, 1:].flatten().astype(float)
        test_values = test_output.values[:, 1:].flatten().astype(float)

        self.assertTrue(np.allclose(example_values, test_values, atol=2.0))

        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_filtering_by_function_list_*.tab"):
            os.remove(fl)

    def test_is_output_correct_for_fishtaco_shapley_value(self):
        """Does FishTaco with no inference produce the correct output for
        the example case with SHAPLEY value calculation?
        """
        print("==============================================================")
        print("Testing compute_contribution_to_DA.py (Shapley value)")
        print("Path to examples:" + FishTacoTestCase.path_to_data)
        print("==============================================================")
        sys.stdout.flush()

        # run fishtaco from shell to make sure it works properly
        subprocess.call('run_fishtaco.py -op fishtaco_out_shapley -mult_hyp '
                        'Bonf' +
                        ' -ta ' +
                        FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab' +
                        ' -fu ' + FishTacoTestCase.path_to_data +
                        '/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab' +
                        ' -l  ' + FishTacoTestCase.path_to_data +
                        '/examples/SAMPLE_vs_CLASS.tab' +
                        ' -gc ' + FishTacoTestCase.path_to_data +
                        '/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab' +
                        ' -assessment multi_taxa -score wilcoxon -na_rep 0 '
                        '-number_of_shapley_orderings_per_taxa 10 -log ' +
                        ' -number_of_permutations 100 '
                        '-single_function_filter K00001 ' +
                        ' -map_function_level none '
                        '-functional_profile_already_corrected_with_musicc',
                        shell=True)

        print("Testing output...")

        # assert that the calculations of the shift contribution of each taxon
        # are the same as the example
        example_output = pd.read_csv(FishTacoTestCase.path_to_data +
                                     '/examples/output/fishtaco_out_shapley_'
                                     'STAT_taxa_contributions_'
                                     'SCORE_wilcoxon_'
                                     'ASSESSMENT_multi_taxa.tab', sep="\t")

        test_output = pd.read_csv('fishtaco_out_shapley_'
                                  'STAT_taxa_contributions_SCORE_'
                                  'wilcoxon_ASSESSMENT_multi_taxa.tab',
                                  sep="\t")

        example_values = example_output.values[:, 1].flatten().astype(float)
        test_values = test_output.values[:, 1].flatten().astype(float)

        self.assertTrue(np.allclose(example_values, test_values, atol=2.0))

        print("Deleting temporary files...")
        for fl in glob.glob("fishtaco_out_shapley_*.tab"):
            os.remove(fl)

###############################################################################

if __name__ == '__main__':
    unittest.main()


