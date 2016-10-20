#!/usr/bin/env python

"""
This is the running script for FishTaco
"""
# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import argparse

# when the test module will be ready:
from fishtaco.compute_contribution_to_DA import main

if __name__ == "__main__":
    # get options from user
    parser = \
        argparse.ArgumentParser(description='Quantify the '
                                            'individual contributions of '
                                            'taxa to shifts observed in '
                                            'functional composition across '
                                            'different sample sets')

    # Required arguments:

    parser.add_argument('-ta', '--taxonomic_abundance_profiles',
                        dest='taxa_abun_file',
                        help='Input file of taxonomic abundance profiles',
                        default=None)

    parser.add_argument('-fu', '--functional_abundance_profiles',
                        dest='function_abun_file',
                        help='Input file of functional abundance profiles',
                        default=None)

    parser.add_argument('-l', '--labels', dest='class_file',
                        help='Input file of label assignment for the two '
                             'sample sets being compared', default=None)

    # Optional arguments:

    parser.add_argument('-gc', '--genomic_content_of_taxa',
                        dest='taxa_to_function_file',
                        help='Input file of genomic content of each taxa',
                        default=None)

    parser.add_argument('-inf', '--perform_inference_of_genomic_content',
                        dest='apply_inference',
                        help='Defines if genome content is inferred (either '
                             'de-novo or prior-based if genomic content is '
                             'also given)', action='store_true')

    parser.add_argument('-label_to_find_enrichment_in', dest='case_label',
                        help='Define sample set label to find enrichment in '
                             '(default: 1)', default='1')

    parser.add_argument('-label_to_find_enrichment_against',
                        dest='control_label',
                        help='Define sample set label to find enrichment '
                             'against (default: 0)', default='0')

    parser.add_argument('-op', '--output_prefix', dest='output_pref',
                        help='Output prefix for result files (default: '
                             'fishtaco_out)', default='fishtaco_out')

    parser.add_argument('-map_function_level', dest='map_function_level',
                        help='Map functions to pathways or modules '
                             '(default: pathway)',
                        default='pathway', choices=['pathway', 'module',
                                                    'none', 'custom'])

    # Advanced arguments:

    parser.add_argument('-map_function_file', dest='map_function_file',
                        help='pathways or modules mapping file (default: use '
                             'internal KEGG file)',
                        default=None)

    parser.add_argument('-perform_inference_on_ko_level',
                        dest='perform_inference_on_ko_level',
                        help='Indicates to perform the inference on the KO '
                             'level (default: use the mapped functional '
                             'level, e.g., pathway)', action='store_true')

    parser.add_argument('-mult_hyp', '--multiple_hypothesis_correction',
                        dest='multiple_hypothesis_correction',
                        help='Multiple hypothesis correction for functional '
                             'enrichment (default: FDR-0.05)',
                        default='FDR-0.05',
                        choices=['Bonf', 'FDR-0.01', 'FDR-0.05', 'FDR-0.1',
                                 'none'])

    parser.add_argument('-max_func', '--maximum_functions_to_analyze',
                        dest='max_da_functions_cases_controls',
                        help='Maximum number of enriched functions to '
                             'consider (default: All)', default=None)

    parser.add_argument('-assessment', '--taxa_assessment_method',
                        dest='taxa_assessment_method',
                        help='The method used when assessing taxa to compute '
                             'individual contributions (default: multi_taxa)',
                        default='multi_taxa',
                        choices=['single_taxa', 'multi_taxa'])

    parser.add_argument('-score', '--score_to_compute',
                        dest='score_to_compute',
                        help='The enrichment score to compute for each '
                             'function (default: wilcoxon)',
                        default='wilcoxon',
                        choices=['t_test', 'mean_diff', 'median_diff',
                                 'wilcoxon', 'log_mean_ratio'])

    parser.add_argument('-max_score', '--max_score_cutoff',
                        dest='max_score_cutoff',
                        help='The maximum score cutoff (for example, '
                             'when dividing by zero) (default: 100)',
                        default='100')

    parser.add_argument('-na_rep', dest='na_rep',
                        help='How to represent NAs in the output (default: '
                             'NA)', default='NA')

    parser.add_argument('-number_of_permutations',
                        dest='number_of_permutations',
                        help='number of permutations (default: 100)',
                        default='100')

    parser.add_argument('-number_of_shapley_orderings_per_taxa',
                        dest='number_of_shapley_orderings_per_taxa',
                        help='number of shapley orderings per taxa '
                             '(default: 5)', default='5')

    # DEPRECATED:
    #parser.add_argument('-use_gc_as_prior',
    # '--use_genomic_content_of_taxa_as_prior', dest='use_t2f_as_prior',
    #                    help='Learn the taxa copy number of each function,
    # using the given genomic content data as prior (default: False)',
    # action='store_true')
    #parser.add_argument('-residual_mode', dest='residual_mode',
    # choices=['as_taxa', 'remove_residual', 'as_baseline'],
    #                    help='How to treat the residual of the functional
    # abundance profile (default: remove_residual)', default='remove_residual')
    #parser.add_argument('-normalization_mode', dest='normalization_mode',
    # choices=['none', 'scale_non_permuted', 'scale_permuted'],
    #                    help='How to normalize the sample after permuting
    # taxa (default: scale_permuted)', default='scale_permuted')
    #parser.add_argument('-permutation_mode', dest='permutation_mode',
    # choices=['independent', 'blocks'],
    #                    help='How to permute the taxa across samples
    # (default: blocks)', default='blocks')

    parser.add_argument('-en', '--enrichment_results', dest='da_result_file',
                        help='Pre-computed functional enrichment results '
                             'from the compute_differential_abundance.py '
                             'script (default: None)', default=None)

    parser.add_argument('-single_function_filter',
                        dest='single_function_filter',
                        help='Limit analysis only to this single function ('
                             'default: None)', default=None)

    parser.add_argument('-multi_function_filter_list',
                        dest='multi_function_filter_list',
                        help='Limit analysis only to these comma-separated '
                             'functions (default: None)', default=None)

    parser.add_argument('-functional_profile_already_corrected_with_musicc',
                        dest='functional_profile_already_corrected_with_musicc',
                        help='Indicates that the functional profile has been '
                             'already corrected with MUSiCC prior to running '
                             'FishTaco (default: False)', action='store_true')

    parser.add_argument('-log', '--log', dest='write_log',
                        help='Write to log file (default: False)',
                        action='store_true')

    given_args = parser.parse_args()
    main(vars(given_args))
