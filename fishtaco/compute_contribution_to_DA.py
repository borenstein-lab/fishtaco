"""
This function quantifies the individual contributions of taxa to the outcome of
functional differential abundance (shift), given the full linear
deconvolution of the functional profile to the underlying taxonomic profile
and the genomic content of each taxon.

Parameters
----------

    args: dictionary
        a dictionary containing the function parameters
        args.['taxa_abun_file']: Input file of taxonomic abundance profiles
        args.['function_abun_file']: Input file of functional abundance
        profiles
        args.['class_file']: Input file of label assignment for the two
        sample sets being compared
        args.['taxa_to_function_file']: Input file of genomic content of
        each taxa
        args.['apply_inference']: Defines if genome content is inferred
        (either de-novo or prior-based if genomic content is also given)
        args.['case_label']: Define sample set label to find enrichment in
        args.['control_label']: Define sample set label to find enrichment
        against
        args.['output_pref']: Output prefix for result files
        args.['map_function_level']: Map functions to pathways or modules
        args.['map_function_file']: pathways or modules mapping file
        args.['perform_inference_on_ko_level']: Indicates to perform the
        inference on the KO level
        args.['multiple_hypothesis_correction']: Multiple hypothesis
        correction for functional enrichment
        args.['max_da_functions_cases_controls']: Maximum number of enriched
        functions to consider
        args.['taxa_assessment_method']: The method used when assessing taxa
        to compute individual contributions
        args.['score_to_compute']: The enrichment score to compute for each
        function
        args.['max_score_cutoff']: The maximum score cutoff (for example,
        when dividing by zero)
        args.['na_rep']: How to represent NAs in the output
        args.['number_of_permutations']: number of permutations
        args.['number_of_shapley_orderings_per_taxa']: number of shapley
        orderings per taxa
        args.['da_result_file']: Pre-computed functional enrichment results
        from the compute_differential_abundance.py script
        args.['single_function_filter']: Limit analysis only to this single
        function
        args.['multi_function_filter_list']: Limit analysis only to these
        comma-separated functions
        args.['functional_profile_already_corrected_with_musicc']: Indicates
        that the functional profile has been already corrected with MUSiCC
        prior to running FishTaco
        args.['write_log']: 'Write to log file

"""

# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function
# imports
import numpy as np
import pandas as pd
import sys
import os
import warnings
import time
import argparse
from math import log
from scipy import stats
from scipy import optimize
from musicc.core import correct_and_normalize
import fishtaco
from fishtaco import compute_pathway_abundance
from fishtaco import compute_differential_abundance
from fishtaco import learn_non_neg_elastic_net_with_prior
from sklearn import cross_validation

__author__ = 'Ohad Manor'
__email__ = 'omanor@gmail.com'
__status__ = "Development"


###############################################################################
# COMPUTE SCORE SUB-FUNCTION
###############################################################################
def compute_score(abundance_cases, abundance_controls, score_to_compute,
                  max_score_cutoff):
    """
    Computes differential abundance score for two abundance vectors

    Parameters
    ----------
    abundance_cases:
        taxa abundance vector for the case samples

    abundance_controls:
        taxa abundance vector for the control samples

    score_to_compute:
        type of differential abundance score to compute

    max_score_cutoff:
        maximum score to return

    Returns
    -------
    score:
        the differential abundance score associted with the given abundance
        vectors

    """

    if score_to_compute == 't_test':
        score, _ = stats.ttest_ind(abundance_cases, abundance_controls)
    elif score_to_compute == 'wilcoxon':
        score, _ = stats.ranksums(abundance_cases, abundance_controls)
    elif score_to_compute == 'mean_diff':
        score = np.mean(abundance_cases) - np.mean(abundance_controls)
    elif score_to_compute == 'median_diff':
        score = np.median(abundance_cases) - np.median(abundance_controls)
    elif score_to_compute == 'log_mean_ratio':
        if np.mean(abundance_controls) == 0:
            score = float(max_score_cutoff)
        elif np.mean(abundance_cases) == 0:
            score = -1 * float(max_score_cutoff)
        elif np.mean(abundance_controls) < 0:
            score = 0
        elif np.mean(abundance_cases) < 0:
            score = 0
        else:
            score = log(np.mean(abundance_cases) / np.mean(abundance_controls))

    return score


###############################################################################
# MAIN FUNCTION
###############################################################################
def main(args):

    # Get the path to the location of the package:
    path_to_data = os.path.dirname(fishtaco.__file__)

    # set some initial settings for the script
    np.set_printoptions(precision=5, suppress=False, linewidth=200)

    ###########################################################################
    # SET VALUES FOR DEPRECATED PARAMETERS
    ###########################################################################
    args['residual_mode'] = 'remove_residual'
    args['normalization_mode'] = 'scale_permuted'
    args['permutation_mode'] = 'blocks'

    print('Given parameters:')
    print(args)

    ###########################################################################
    # CHECK FOR PARAMETER CONTRADICTIONS
    ###########################################################################
    if args['residual_mode'] == 'as_taxa' and \
                    args['normalization_mode'] != 'none':
        sys.exit('Error: residual_mode==as_taxa cannot be used with '
                 'normalization_mode!=none')

    if args['permutation_mode'] == 'independent' and \
                    args['normalization_mode'] != 'none':
        sys.exit('Error: permutation_mode==independent cannot be used with '
                 'normalization_mode!=none')

    if args['normalization_mode'] == 'scale_non_permuted' and \
                    args['taxa_assessment_method'] != 'multi_taxa':
        sys.exit('Error: currently normalization_mode==scale_non_permuted '
                 'can only be used with taxa_assessment_method==multi_taxa')

    if args['normalization_mode'] == 'scale_permuted' and \
                    args['taxa_assessment_method'] != 'multi_taxa' and \
                    args['taxa_assessment_method'] != 'single_taxa':
        sys.exit('Error: currently normalization_mode==scale_permuted can '
                 'only be used with taxa_assessment_method==multi_taxa or '
                 'taxa_assessment_method==single_taxa')

    ###########################################################################
    # CREATE OUTPUT SUFFIX
    ###########################################################################
    output_suffix = '_SCORE_' + args['score_to_compute'] + '_ASSESSMENT_' + \
                    args['taxa_assessment_method'] + '.tab'

    ###########################################################################
    # OPEN LOG FILE IF REQUESTED
    ###########################################################################
    if 'write_log' in args.keys() and args['write_log']:
        with open(args['output_pref'] + '_STAT_run_log' +
                          output_suffix, 'w') as f:
            f.write("# " + sys.argv[0] + " " + str(args) + '\n')

    ###########################################################################
    # INPUT
    ###########################################################################
    print("Reading input files...")

    # read taxonomic abundance file
    if args['taxa_abun_file'] is not None:
        if not os.path.isfile(args['taxa_abun_file']):
            sys.exit('Error: Input file "' + args['taxa_abun_file'] +
                     '" does not exist')
        original_taxa_abun_data = pd.read_table(args['taxa_abun_file'],
                                                index_col=0, dtype={0: str})

        if np.sum(np.isnan(original_taxa_abun_data.values)) > 0:
            sys.exit('Error: Taxa abundance contains NaN')

        if args['normalization_mode'] == 'scale_non_permuted' or \
                        args['normalization_mode'] == 'scale_permuted':
            # normalize taxa abundance of all taxa to 1
            normalized_taxa_abundance = original_taxa_abun_data.values / \
                                        np.sum(original_taxa_abun_data.values,
                                               axis=0)
            original_taxa_abun_data = \
                pd.DataFrame(data=normalized_taxa_abundance,
                             index=original_taxa_abun_data.index.values,
                             columns=original_taxa_abun_data.columns.values)

    else:
        sys.exit('Error: No input taxa abundance given to script')

    # read function abundance file
    if args['function_abun_file'] is not None:
        if not os.path.isfile(args['function_abun_file']):
            sys.exit('Error: Input file "' + args['function_abun_file'] +
                     '" does not exist')

        if 'functional_profile_already_corrected_with_musicc' in args.keys() \
                and args['functional_profile_already_corrected_with_musicc']:

            function_abun_data = pd.read_table(args['function_abun_file'],
                                               index_col=0)

        else:  # if needed, correct the functional profile using MUSiCC

            args_for_running_musicc = \
                {'input_file': args['function_abun_file'],
                 'output_file': args['output_pref'] +
                                '_STAT_function_abundance_MUSiCC_corrected' +
                                output_suffix,
                 'input_format': 'tab', 'output_format': 'tab',
                 'musicc_inter': True, 'musicc_intra': 'learn_model',
                 'compute_scores': False, 'verbose': False}

            try:
                correct_and_normalize(args_for_running_musicc)

            except:
                print("Unexpected error while running MUSiCC on functional "
                      "data (make sure you have KEGG orthology groups (KOs) "
                      "as rows):")
                print(sys.exc_info())
                print("-----------------------------------------------------")
                sys.exit()

            # read the corrected data into the variable
            function_abun_data = \
                pd.read_table(args['output_pref'] +
                              '_STAT_function_abundance_MUSiCC_corrected' +
                              output_suffix, index_col=0)

        # save original function abundance if needed later
        original_function_abun_data = function_abun_data

        if np.sum(np.isnan(function_abun_data.values)) > 0:
            sys.exit('Error: Function abundance contains NaN')
    else:
        print('No input of functional abundance given to FishTaco, '
              'predicting from taxonomic abundance and genomic content...')

    # read taxa_to_function (functional genomic content of each taxon) file
    if args['taxa_to_function_file'] is not None:
        if not os.path.isfile(args['taxa_to_function_file']):
            sys.exit('Error: Input file "' + args['taxa_to_function_file'] +
                     '" does not exist')
        taxa_to_function_data = pd.read_table(args['taxa_to_function_file'],
                                              index_col=0, dtype={0: str})
        taxa_to_function_data.index.name = "Taxa"

        # save original taxa to function if needed later
        original_taxa_to_function_data = taxa_to_function_data
        original_taxa_to_function_data.index.name = "Taxa"

        if np.sum(np.isnan(taxa_to_function_data.values)) > 0:
            sys.exit('Error: Taxa to function data contains NaN')

        # if we did no receive as input functional profiles,
        # we can predict them from the taxonomic abundance and the genomic
        # content:
        if args['function_abun_file'] is None:
            args_for_predicting_function_from_taxa_and_content = \
                {'ko_abun_pd': original_taxa_abun_data,
                 'ko_to_pathway_pd': original_taxa_to_function_data,
                 'output_pd': pd.DataFrame(),
                 'mapping_method': 'naive', 'compute_method': 'sum',
                 'transpose_ko_abundance': False, 'transpose_output': False,
                 'verbose': False}

            compute_pathway_abundance.main(
                args_for_predicting_function_from_taxa_and_content)

            function_abun_data = \
                args_for_predicting_function_from_taxa_and_content['output_pd']

            # save original function abundance if needed later
            original_function_abun_data = function_abun_data

    else:
        if args['function_abun_file'] is None:
            sys.exit('Error: No input of functional abundance or genomic '
                     'content given to FishTaco, exiting...')
        # with no given genomic content, we need to perform inference:
        elif 'apply_inference' in args.keys() and args['apply_inference']:
            print('No input of genomic content given to FishTaco, inferring '
                  'the mapping of taxa to functions from taxonomic and '
                  'functional profiles')
        else:
            sys.exit('Error: No input taxa to function data given to script,'
                     'and inference not requested')

    # check if we have to map the functions to a higher level
    if args['map_function_level'] != 'none':

        # read mapping file from function to higher level
        if args['map_function_file'] is not None:
            if not os.path.isfile(args['map_function_file']):
                sys.exit('Error: Input file "' + args['map_function_file'] +
                         '" does not exist')

            function_mapping_file = args['map_function_file']

        else:
            if args['map_function_level'] == 'pathway':
                function_mapping_file = \
                    path_to_data + \
                    '/data/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15.tab'

            elif args['map_function_level'] == 'module':
                function_mapping_file = \
                    path_to_data + \
                    '/data/KOvsMODULE_BACTERIAL_KEGG_2013_07_15.tab'

            else:  # custom
                sys.exit('Error: No custom mapping file given')

        # map functions to pathway level
        print("Mapping functions to pathway/module level...")
        function_to_pathway_mapping = \
            pd.read_table(function_mapping_file, index_col=0, dtype={0: str})

        # filter out functions from the mapping file that we don't have in our
        # functional profile:
        functions_in_mapping_with_abundance = \
            np.intersect1d(function_to_pathway_mapping.index.values,
                           original_function_abun_data.index.values)
        if len(functions_in_mapping_with_abundance) == 0:
            sys.exit('Error: None of the functions in the mapping file appear '
                     'in the function abundance file')

        function_to_pathway_mapping = \
            function_to_pathway_mapping.loc[functions_in_mapping_with_abundance]

        args_for_mapping_function_to_higher = {'ko_abun_pd': function_abun_data,
                                               'ko_to_pathway_pd': function_to_pathway_mapping,
                                               'output_pd': pd.DataFrame(),
                                               'mapping_method': 'naive',
                                               'compute_method': 'sum',
                                               'transpose_ko_abundance': False,
                                               'transpose_output': False, 'verbose': False}

        compute_pathway_abundance.main(args_for_mapping_function_to_higher)

        function_abun_data = args_for_mapping_function_to_higher['output_pd']

        if np.sum(np.isnan(function_abun_data.values)) > 0:
            sys.exit('Error: Function abundance after mapping contains NaN')

        # map the taxa to functions file to the pathway level
        if args['taxa_to_function_file'] is not None:
            args_for_mapping_function_to_higher = {'ko_abun_pd': original_taxa_to_function_data,
                                                   'ko_to_pathway_pd': function_to_pathway_mapping,
                                                   'output_pd': pd.DataFrame(),
                                                   'mapping_method': 'naive',
                                                   'compute_method': 'sum',
                                                   'transpose_ko_abundance': True,
                                                   'transpose_output': True,
                                                   'verbose': False}

            compute_pathway_abundance.main(args_for_mapping_function_to_higher)

            taxa_to_function_data = args_for_mapping_function_to_higher['output_pd']

        print("Done.")

    # read class labels file
    if args['class_file'] is not None:
        if not os.path.isfile(args['class_file']):
            sys.exit('Error: Input file "' + args['class_file'] + '" does not exist')
        class_data = pd.read_table(args['class_file'], index_col=0, dtype=str)

    else:
        sys.exit('Error: No input class data given to script')

    print("Done.")

    # first, reduce all input files to the same samples, sorted
    print("Reducing taxa, function, and class data to contain the exact "
          "same set of samples...")
    samples = np.sort(np.intersect1d(np.intersect1d(function_abun_data.columns.values,
                                                    original_taxa_abun_data.columns.values),
                                     class_data.index.values))
    function_abun_data = function_abun_data[samples]
    original_function_abun_data = original_function_abun_data[samples]
    original_taxa_abun_data = original_taxa_abun_data[samples]
    class_data = class_data.loc[samples]
    #print(class_data)
    print("Done.")

    # If needed, filter our functions to the single function or multi functions given
    # note that we pass as a list to ensure we get a DataFrame back!
    if args['single_function_filter'] is not None:
        print("Filtering for a single function: " + args['single_function_filter'])
        function_abun_data = function_abun_data.loc[[args['single_function_filter']], :]
        print("Done.")

    # note that we pass as a list to ensure we get a DataFrame back!
    if args['multi_function_filter_list'] is not None:
        print("Filtering for the following functions: " +
              args['multi_function_filter_list'])
        multi_function_filter_as_array = args['multi_function_filter_list'].split(',')
        function_abun_data = function_abun_data.loc[multi_function_filter_as_array, :]
        print("Done.")

    ###########################################################################
    # TEST THAT WE HAVE SUFFICIENT SAMPLES (AT LEAST 3 FROM EACH CLASS)
    ###########################################################################
    # define controls and cases
    number_of_samples = function_abun_data.shape[1]
    controls = (class_data.values.reshape(number_of_samples) == args['control_label'])
    cases = (class_data.values.reshape(number_of_samples) == args['case_label'])
    print("#cases = " + str(sum(controls)) + ", #controls = " + str(sum(cases)))

    if sum(controls) < 3 or sum(cases) < 3:
        print("Error: Cases or Controls have less than 3 samples, exiting...")
        exit()

    ###########################################################################
    # TAXA DIFFERENTIAL ABUNDANCE (SHIFT)
    ###########################################################################
    # compute a differential abundance score for each taxa:
    if args['score_to_compute'] == 't_test' or args['score_to_compute'] == \
            'mean_diff' or args['score_to_compute'] == 'log_mean_ratio':
        da_method = 'ttest'
    else:  # score_to_compute == median_diff or score_to_compute == wilcoxon
        da_method = 'wilcoxon'

    print("Computing a differential abundance score for each taxa...")
    args_for_taxa_differential_abundance = {'input_pd': original_taxa_abun_data,
                                            'class_pd': class_data,
                                            'output_pd': pd.DataFrame(),
                                            'method': da_method,
                                            'control_label': args['control_label'],
                                            'case_label': args['case_label'],
                                            'class_header': True,
                                            'verbose': False}

    compute_differential_abundance.main(args_for_taxa_differential_abundance)
    #print(args_for_TAXA_differential_abundance['output_pd'])
    taxa_diff_abun_scores = args_for_taxa_differential_abundance['output_pd']
    print("Done.")

    ###########################################################################
    # FUNCTION DIFFERENTIAL ABUNDANCE (SHIFT)
    ###########################################################################
    # If we are given a pre-computed functional shift file, read the p-values
    # and filter using the threshold:
    if 'da_result_file' in args.keys() and args['da_result_file'] is not None:
        print("Using function differential abundance scores given...")
        da_scores = pd.read_table(args['da_result_file'], index_col=0)
        if args['multiple_hypothesis_correction'] == 'none':
            da_functions = da_scores.index.values
        else:
            da_functions = da_scores[da_scores[args['multiple_hypothesis_correction']] > 0].index.values

    else:  # run the functional shift analysis on the functional profiles we got:
        print("Computing a differential abundance score for each fucntion...")
        args_for_differential_abundance = {'input_pd': function_abun_data,
                                           'class_pd': class_data,
                                           'output_pd': pd.DataFrame(),
                                           'method': da_method,
                                           'control_label': args['control_label'],
                                           'case_label': args['case_label'],
                                           'class_header': True,
                                           'verbose': False}

        compute_differential_abundance.main(args_for_differential_abundance)
        da_scores = args_for_differential_abundance['output_pd']
        if args['multiple_hypothesis_correction'] == 'none':
            da_functions = da_scores.index.values
        else:
            da_functions = da_scores[da_scores[args['multiple_hypothesis_correction']] > 0].index.values

    print("Done.")
    # for testing:
    # print(da_scores)
    # print(da_functions)

    ###########################################################################
    # SELECT ONLY FUNCTIONS THAT ARE ENRICHED IN THE SPECIFIC SET OF SAMPLES
    # WITH THE CASES LABEL
    ###########################################################################
    print("Selecting only functions that are enriched in samples with the label: " +
          args['case_label'] + "...")
    da_scores.loc[np.isnan(da_scores['singLogP']), 'singLogP'] = 0
    enriched_in_cases = da_scores.index.values[np.where(da_scores['singLogP'] > 0)]
    da_functions = np.intersect1d(enriched_in_cases, da_functions)
    print("Done.")

    ###########################################################################
    # SELECT ONLY TOP DA FEATURES FOR FURTHER ANALYSIS
    ###########################################################################
    if 'max_da_functions_cases_controls' in args.keys() and \
                    args['max_da_functions_cases_controls'] is not None:
        print("Selecting only " + args['max_da_functions_cases_controls'] +
              " differentially abundant functions...")
        da_scores.loc[np.isnan(da_scores['singLogP']), 'singLogP'] = 0
        controls_top_functions = \
            da_scores.sort('singLogP')[0:int(args['max_da_functions_cases_controls'])].index.values
        cases_top_functions = \
            da_scores.sort('singLogP')[da_scores.shape[0] -
                                       int(args['max_da_functions_cases_controls']):da_scores.shape[0]].index.values
        da_functions = np.intersect1d(np.union1d(cases_top_functions,
                                                 controls_top_functions), da_functions)
        print("Done.")

    ###########################################################################
    # REDUCE TO DA FUNCTIONS
    ###########################################################################
    # now, reduce the input functions to be only the DA ones and sort them
    if args['taxa_to_function_file'] is not None:
        functions = np.sort(np.intersect1d(np.intersect1d(function_abun_data.index.values,
                                                          taxa_to_function_data.columns.values),
                                           da_functions))
        function_abun_data = function_abun_data.loc[functions]
        functions_da_scores = da_scores.loc[functions]
        taxa_to_function_data = taxa_to_function_data[functions]
    else:
        functions = np.sort(np.intersect1d(function_abun_data.index.values, da_functions))
        function_abun_data = function_abun_data.loc[functions]
        functions_da_scores = da_scores.loc[functions]

    # now, reduce the input taxa to be the same in the two input files and sort
    # them
    if args['taxa_to_function_file'] is not None:
        taxa = np.sort(np.intersect1d(original_taxa_abun_data.index.values.astype(str),
                                      taxa_to_function_data.index.values.astype(str)))
        original_taxa_abun_data = original_taxa_abun_data.loc[taxa]
        taxa_to_function_data = taxa_to_function_data.loc[taxa]
        original_taxa_to_function_data = original_taxa_to_function_data.loc[taxa]
        taxa_diff_abun_scores = taxa_diff_abun_scores.loc[taxa]
    else:
        taxa = np.sort(original_taxa_abun_data.index.values)
        original_taxa_abun_data = original_taxa_abun_data.loc[taxa]
        taxa_diff_abun_scores = taxa_diff_abun_scores.loc[taxa]

    num_of_da_functions = function_abun_data.shape[0]
    num_of_taxa = original_taxa_abun_data.shape[0]

    print("#DA functions:" + str(num_of_da_functions) + " #Taxa:" +
          str(num_of_taxa) + " #Samples:" + str(number_of_samples))
    if 'write_log' in args.keys() and args['write_log']:
        with open(args['output_pref'] + '_STAT_run_log' + output_suffix, 'a') as f:
            f.write("#DA functions:" + str(num_of_da_functions) + " #Taxa:" +
                    str(num_of_taxa) + " #Samples:" + str(number_of_samples) + "\n")

    ###########################################################################
    # IF NEEDED/REQUESTED, LEARN A PREDICTED COPY-NUMBER FOR EACH SPECIES FOR
    # EACH DA FUNCTION for each DA function, learn a non-negative
    # elastic-net model from taxa, using the genomic content as prior if given.
    # NOTE: We are using the elastic net only as a feature selection step,
    # and use NNLS to fit the exact copy-number
    ###########################################################################
    if num_of_da_functions > 0 and 'apply_inference' in args.keys() and \
            args['apply_inference']:

        print("Inferring the genomic content of each taxa...")

        # if we need to perform the inference on KO level, we need to identify
        # KOs that are mapped to DA pathways, and perform the inference only
        #  on them
        if args['map_function_level'] != 'none' and \
                        'perform_inference_on_ko_level' in args.keys() and \
                args['perform_inference_on_ko_level']:

            func_map_sum = np.sum(function_to_pathway_mapping.loc[:, functions].values, axis=1)

            functions_to_infer = function_to_pathway_mapping.index.values[func_map_sum > 0]

            num_of_functions_to_infer = len(functions_to_infer)

        # we are inferring content for higher level functions (e.g., pathways),
        # or there is no mapping done at all
        else:
            num_of_functions_to_infer = num_of_da_functions
            functions_to_infer = functions

        num_cv = 5

        all_functions_mean_cv_test_rsqr = np.zeros(num_of_functions_to_infer)
        all_functions_global_cv_test_stats = np.zeros((num_of_functions_to_infer, 3))
        all_functions_mean_cv_validation_rsqr = np.zeros(num_of_functions_to_infer)
        all_functions_prediction_on_test_uncentered = np.zeros((number_of_samples,
                                                                num_of_functions_to_infer))
        all_taxa_to_function_weights = np.zeros((num_of_taxa, num_of_functions_to_infer))

        for i in range(num_of_functions_to_infer):

            curr_function = functions_to_infer[i]
            print(i, ":", curr_function)

            # if we are inferring from KO level, then use the original
            # functional abundance data
            if 'perform_inference_on_ko_level' in args.keys() and \
                    args['perform_inference_on_ko_level']:
                function_abun_for_curr_func = \
                    original_function_abun_data.loc[curr_function].values

            # we have a mapping to a higher functional level and we are
            # performing inference for that level
            else:
                function_abun_for_curr_func = \
                    function_abun_data.loc[curr_function].values

            function_abun_non_zero = np.sum(function_abun_for_curr_func != 0)

            # create cross-validation indices
            number_of_cases = np.sum((class_data.values.reshape(number_of_samples) == args['case_label']))
            number_of_controls = np.sum((class_data.values.reshape(number_of_samples) == args['control_label']))
            cases_indices = (class_data.values.reshape(number_of_samples) == args['case_label'])
            control_indices = (class_data.values.reshape(number_of_samples) == args['control_label'])
            k_fold_cases = cross_validation.KFold(number_of_cases, n_folds=num_cv, shuffle=True)
            k_fold_controls = cross_validation.KFold(number_of_controls, n_folds=num_cv, shuffle=True)
            merged_test_indices = np.zeros((number_of_samples, 5))
            for k, (train_cases, test_cases) in enumerate(k_fold_cases):
                cases_integer_indices = np.nonzero(cases_indices)[0]
                merged_test_indices[cases_integer_indices[test_cases], k] = 1
            for k, (train_control, test_control) in enumerate(k_fold_controls):
                control_integer_indices = np.nonzero(control_indices)[0]
                merged_test_indices[control_integer_indices[test_control], k] = 1

            # If we don't have a genomic content file, do de novo inference
            if args['taxa_to_function_file'] is None:
                params = {'covariates_prior': np.ones(num_of_taxa)}

            else:  # if we were given a genomic content file, we use it as a prior
                if 'perform_inference_on_ko_level' in args.keys() and \
                        args['perform_inference_on_ko_level']:
                    function_index_in_taxa_to_function = \
                        original_taxa_to_function_data.columns.get_loc(curr_function)
                    params = {'covariates_prior':
                                  original_taxa_to_function_data.values[:, function_index_in_taxa_to_function]}

                else:
                    function_index_in_taxa_to_function = \
                        taxa_to_function_data.columns.get_loc(curr_function)
                    params = {'covariates_prior':
                                  taxa_to_function_data.values[:, function_index_in_taxa_to_function]}

                # set the minimal prior to 10% of minimal value instead of 0
                if np.sum(params['covariates_prior'] > 0) > 0:
                    minimal_prior_value = np.min(params['covariates_prior'][params['covariates_prior'] > 0]) * 0.1
                else:
                    minimal_prior_value = 0.1
                params['covariates_prior'][params['covariates_prior'] == 0] = minimal_prior_value

            if np.max(params['covariates_prior']) > 0 and function_abun_non_zero > 1:

                test_rsqr = np.zeros(num_cv)
                mean_validation_rsqr = np.zeros(num_cv)
                model_weights = np.zeros((num_of_taxa, num_cv))
                prediction_on_test = np.zeros(number_of_samples)

                for c in range(num_cv):

                    test_indices = np.nonzero(merged_test_indices[:, c])[0]
                    train_indices = np.nonzero(merged_test_indices[:, c] == 0)[0]

                    cov_train = original_taxa_abun_data.transpose().values[train_indices, :]

                    cov_test = original_taxa_abun_data.transpose().values[test_indices, :]

                    res_train = function_abun_for_curr_func[train_indices]

                    res_test = function_abun_for_curr_func[test_indices]

                    # with NO class subfeatures
                    params['class_subfeatures'] = None
                    if np.sum(np.isnan(cov_train)) > 0:
                        print(curr_function)
                        print(cov_train)
                        sys.exit('Error: cov_train contains NaN')
                    if np.sum(np.isnan(res_train)) > 0:
                        print(curr_function)
                        print(res_train)
                        sys.exit('Error: res_train contains NaN')

                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", category=UserWarning)
                        enet, validation_rsqr = learn_non_neg_elastic_net_with_prior.learn(cov_train, res_train, params)

                    mean_validation_rsqr[c] = np.mean(validation_rsqr)
                    model_weights[:, c] = enet.coef_
                    prediction_on_test[test_indices], test_rsqr[c] = \
                        learn_non_neg_elastic_net_with_prior.test(enet, cov_test, res_test, params)

                    # now use NNLS to learn the final weights, using only
                    # features that got non-zero weights in the model

                    non_zero_features = enet.coef_ > 0

                    # FIT A NNLS MODEL TO THE NON-ZERO FEATURES FROM THE
                    # ELASTIC NET
                    if sum(non_zero_features) > 0:
                        try:
                            nnls_weights, _ = \
                                optimize.nnls(original_taxa_abun_data.values[non_zero_features, :].T,
                                              function_abun_for_curr_func)

                        except ValueError:
                            print("Oops, an error occurred during inference!")
                            print("Number of non-zero features was: " +
                                  str(sum(non_zero_features)))
                            print("Feature values were: " +
                                  str(original_taxa_abun_data.values[non_zero_features, :].T))
                            print("Response was: " +
                                  str(function_abun_for_curr_func))
                            sys.exit()

                        model_weights[non_zero_features, c] = nnls_weights

                    prediction_on_test[test_indices] = np.dot(original_taxa_abun_data.values[:, test_indices].T, model_weights[:, c])
                    all_functions_prediction_on_test_uncentered[test_indices, i] = prediction_on_test[test_indices]

                all_taxa_to_function_weights[:, i] = np.median(model_weights, axis=1)
                all_functions_mean_cv_test_rsqr[i] = np.mean(test_rsqr)

                all_functions_global_cv_test_stats[i, 0] = 1 - (np.sum((prediction_on_test - function_abun_for_curr_func) ** 2) /
                                                                np.sum((function_abun_for_curr_func - np.mean(function_abun_for_curr_func)) ** 2))
                all_functions_global_cv_test_stats[i, 1], _ = stats.pearsonr(prediction_on_test, function_abun_for_curr_func)
                all_functions_global_cv_test_stats[i, 2], _ = stats.spearmanr(prediction_on_test, function_abun_for_curr_func)
                all_functions_mean_cv_validation_rsqr[i] = np.mean(mean_validation_rsqr)

            else:
                all_taxa_to_function_weights[:, i] = 0
                all_functions_mean_cv_test_rsqr[i] = 0
                all_functions_global_cv_test_stats[i, 0] = 0
                all_functions_global_cv_test_stats[i, 1] = 0
                all_functions_global_cv_test_stats[i, 2] = 0
                all_functions_mean_cv_validation_rsqr[i] = 0

        # now, replace the original taxa-to-function data by the learned weights and continue
        if args['map_function_level'] != 'none' and \
                        'perform_inference_on_ko_level' in args.keys() and \
                args['perform_inference_on_ko_level']:

            original_taxa_to_function_data = pd.DataFrame(data=all_taxa_to_function_weights,
                                                          index=original_taxa_abun_data.index.values,
                                                          columns=functions_to_infer)
            original_taxa_to_function_data.index.name = 'Taxa'

            args_for_mapping_function_to_higher = {'ko_abun_pd': original_taxa_to_function_data,
                                                   'ko_to_pathway_pd': function_to_pathway_mapping,
                                                   'output_pd': pd.DataFrame(),
                                                   'mapping_method': 'naive',
                                                   'compute_method': 'sum',
                                                   'transpose_ko_abundance': True,
                                                   'transpose_output': True,
                                                   'verbose': False}

            compute_pathway_abundance.main(args_for_mapping_function_to_higher)

            taxa_to_function_data = args_for_mapping_function_to_higher['output_pd']
            taxa_to_function_data = taxa_to_function_data[functions]

        else:

            original_taxa_to_function_data = pd.DataFrame(data=all_taxa_to_function_weights,
                                                          index=original_taxa_abun_data.index.values,
                                                          columns=functions_to_infer)
            original_taxa_to_function_data.index.name = 'Taxa'

            taxa_to_function_data = pd.DataFrame(data=all_taxa_to_function_weights,
                                                 index=original_taxa_abun_data.index.values,
                                                 columns=functions_to_infer)
            taxa_to_function_data.index.name = 'Taxa'

    ###########################################################################
    # FIRST, GIVEN THE TAXA-TO-FUNCTION COPY-NUMBER/WEIGHTS, ADD A NEW
    # "TAXA"  TO THE MATRICES WITH THE NAME "UNKNOWN". THIS TAXA WILL BEHAVE
    # EXACTLY LIKE ALL THE REST FOR COMPUTATION OF CONTRIBUTION, EXCEPT
    # FROM  THE FACT THAT IT WILL HAVE A DIFFERENT ABUNDANCE FOR EVERY
    # FUNCTION, SINCE THE RESIDUAL OF EVERY FUNCTION IS DIFFERENT.
    # THEREFORE, WE ARE NOT REMOVING THE RESIDUAL, BUT RATHER NAMING IT
    # "UNKNOWN", AND CALCULATING EVERYTHING IN THE FULL CONTEXT OF ALL THE
    # TAXA AND THE RESIDUAL NOTE THAT ALSO WITH INFERENCE, WE USE THE FINAL
    # AND FULL MODEL (I.E., THE MEDIAN OF THE LEARNED WEIGHTS) TO PREDICT
    # THE RESIDUAL, SINCE WE WANT THE THE RESIDUAL TO DESCRIBE THE ACTUAL
    # ERROR WHEN WE MULTIPLY THE TAXA ABUNDANCE BY THE INFERRED COPY NUMBERS
    ###########################################################################
    if num_of_da_functions > 0:
        if 'apply_inference' in args.keys() and args['apply_inference']:
            predicted_function_abundance = np.dot(original_taxa_abun_data.values.T,
                                                  taxa_to_function_data.values)
        else:
            predicted_function_abundance = np.dot(original_taxa_abun_data.values.T,
                                                  taxa_to_function_data.values)

    else:
        predicted_function_abundance = function_abun_data.values.T

    # create the residual matrix of functions vs samples for all functions
    residual_function_vs_sample = (function_abun_data.values.T -
                                   predicted_function_abundance).T

    ###########################################################################
    # FOR EACH FUNCTION, COMPUTE THE AGREEMENT BETWEEN THE
    # PREDICTED_FUNCTION ABUNDANCE AND THE REAL FUNCTION ABUNDANCE
    # WE COMPUTE 7 VALUES:
    # 1) THE R^2
    # 2+3) THE PEARSON CORRELATION + P-VALUE
    # 4+5) THE SPEARMAN CORRELATION + P-VALUE
    # 6+7)THE MEAN AND STD OF THE ABSOLUTE VALUE DIFFERENCE
    ###########################################################################
    print("Calculating agreement between metagenome-based and taxa-based "
          "functional profiles...")
    predicted_function_agreement = np.zeros((num_of_da_functions, 7))

    for i in range(num_of_da_functions):

        sos_residual = np.sum((predicted_function_abundance[:, i] -
                               function_abun_data.values[i, :]) ** 2)

        sos_original = np.sum((function_abun_data.values[i, :] -
                               np.mean(function_abun_data.values[i, :])) ** 2)

        predicted_function_agreement[i, 0] = 1 - (sos_residual / sos_original)

        if np.sum(predicted_function_abundance[:, i]) == 0:
            predicted_function_agreement[i, 1] = 0
            predicted_function_agreement[i, 2] = 1
            predicted_function_agreement[i, 3] = 0
            predicted_function_agreement[i, 4] = 1
        else:
            predicted_function_agreement[i, 1], predicted_function_agreement[i, 2] = \
                stats.pearsonr(predicted_function_abundance[:, i], function_abun_data.values[i, :])
            predicted_function_agreement[i, 3], predicted_function_agreement[i, 4] = \
                stats.spearmanr(predicted_function_abundance[:, i], function_abun_data.values[i, :])

        predicted_function_agreement[i, 5] = np.mean(abs(predicted_function_abundance[:, i] -
                                                         function_abun_data.values[i, :]))

        predicted_function_agreement[i, 6] = np.std(abs(predicted_function_abundance[:, i] -
                                                        function_abun_data.values[i, :]))


    ###########################################################################
    # FOR EACH FUNCTION, COMPUTE THE PREDICTED DA GIVEN THE PREDICTED FUNCTION
    # ABUNDANCE
    ###########################################################################
    print("Calculating differential abundance values for each function...")
    original_stat_value = np.zeros(num_of_da_functions)
    predicted_da_stat_value = np.zeros(num_of_da_functions)

    for i in range(num_of_da_functions):
        original_stat_value[i] = compute_score(function_abun_data.values[i, cases],
                                               function_abun_data.values[i, controls],
                                               args['score_to_compute'],
                                               args['max_score_cutoff'])
        predicted_da_stat_value[i] = compute_score(predicted_function_abundance[cases, i],
                                                   predicted_function_abundance[controls, i],
                                                   args['score_to_compute'],
                                                   args['max_score_cutoff'])

    ###########################################################################
    # IF THE USER SELECTED TO REMOVE THE RESIDUAL, THEN SET THE FUNCTIONAL
    # PROFILE TO BE THE PREDICTED ONE (I.E., THE PSEUDO-METAGENOME,
    # AND RESET THE RESIDUAL TO BE ALL ZEROS (AND THUS IT WILL HAVE NO EFFECT)
    ###########################################################################
    if 'residual_mode' in args.keys() and args['residual_mode'] == 'remove_residual':
        function_abun_data = pd.DataFrame(data=predicted_function_abundance.T,
                                          index=function_abun_data.index.values,
                                          columns=function_abun_data.columns.values)
        residual_function_vs_sample = np.zeros((num_of_da_functions, number_of_samples))

    ###########################################################################
    # CREATE PERMUTATION MATRICES FOR RESIDUAL AND REAL TAXA
    ###########################################################################
    number_of_permutations = int(args['number_of_permutations'])
    no_residual_permuted_taxa_abundance_matrices = np.zeros((number_of_permutations,
                                                             number_of_samples, num_of_taxa))
    permuted_residual_matrices_function_vs_samples = np.zeros((number_of_permutations,
                                                               num_of_da_functions, number_of_samples))

    print("Creating " + str(number_of_permutations) + " permutations...")

    for p in range(number_of_permutations):
        # create the permuted matrix for the taxa abundance
        if args['permutation_mode'] == 'independent':
            for j in range(num_of_taxa):
                no_residual_permuted_taxa_abundance_matrices[p, :, j] = \
                    np.random.permutation(original_taxa_abun_data.values[j, :])

        # permute in blocks, so simply permute the rows keeping all taxa
        # together in each sample
        else:
            sample_permutation_index = np.random.permutation(number_of_samples)
            no_residual_permuted_taxa_abundance_matrices[p] = \
                original_taxa_abun_data.values[:, sample_permutation_index].T

        # if needed, create the permuted matrix for the residuals of all
        # functions.
        # note that if we use as_baseline, we never permute the residual,
        # so we just create copies of the residual as the permuted matrices
        if 'residual_mode' in args.keys() and args['residual_mode'] != 'remove_residual':
            for i in range(num_of_da_functions):
                if 'residual_mode' in args.keys() and args['residual_mode'] == 'as_taxa':
                    permuted_residual_matrices_function_vs_samples[p, i, :] = \
                        np.random.permutation(residual_function_vs_sample[i, :])
                else:  # as_baseline
                    if ['taxa_assessment_method'] == 'multi_taxa':
                        permuted_residual_matrices_function_vs_samples[p, i, :] = \
                            residual_function_vs_sample[i, :]
                    else:
                        permuted_residual_matrices_function_vs_samples[p, i, :] = \
                            np.random.permutation(residual_function_vs_sample[i, :])

    print("Done.")

    ###########################################################################
    # IF WE ARE NOT REMOVING THE RESIDUAL, THEN UPDATE NUM_OF_TAXA TO INCLUDE
    # UNKNOWN
    ###########################################################################
    if 'residual_mode' in args.keys() and args['residual_mode'] != 'remove_residual':
        num_of_taxa += 1

    ###########################################################################
    # IF WE ARE USING SHAPLEY, CREATE RANDOM ORDERINGS
    ###########################################################################
    if args['taxa_assessment_method'] == 'multi_taxa':

        number_of_orderings = \
            num_of_taxa * int(args['number_of_shapley_orderings_per_taxa'])

        # create the random orderings to be used by all the permuted matrices
        random_taxa_orderings = np.zeros((number_of_orderings, num_of_taxa))
        for o in range(number_of_orderings):
            random_taxa_orderings[o, :] = np.random.permutation(num_of_taxa)

        random_taxa_orderings = random_taxa_orderings.astype(int)

        # write random orderings to output file
        orderings_pd = pd.DataFrame(data=random_taxa_orderings,
                                    index=(np.arange(number_of_orderings)+1),
                                    columns=(np.arange(num_of_taxa)+1))
        orderings_pd.index.name = 'Ordering'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            orderings_pd.to_csv(args['output_pref'] +
                                '_STAT_shapley_orderings' + output_suffix,
                                sep='\t', na_rep=args['na_rep'])

    ###########################################################################
    # ANALYSIS USING THE STATISTIC, THE GIVEN OR INFERRED COPY-NUMBER OF TAXA
    # IN FUNCTIONS, WITH PERMUTATION OF TAXA
    ###########################################################################
    contribution_matrix = np.zeros((num_of_taxa, num_of_da_functions))
    mean_stat_value = np.zeros((num_of_taxa, num_of_da_functions))
    median_stat_value = np.zeros((num_of_taxa, num_of_da_functions))
    std_stat_value = np.zeros((num_of_taxa, num_of_da_functions))

    for i in range(num_of_da_functions):  # num_of_da_functions
        start = time.time()

        # create the weights vector for this function, and add "1" for the
        # "unknown"
        if 'residual_mode' in args.keys() and \
                        args['residual_mode'] != 'remove_residual':
            weights_of_this_function = \
                np.hstack((taxa_to_function_data.values[:, i], 1))
        else:
            weights_of_this_function = taxa_to_function_data.values[:, i]

        # create the taxa abundance data for this function, including the
        # residual as the "unknown" taxa (if not removed)
        if 'residual_mode' in args.keys() and \
                        args['residual_mode'] != 'remove_residual':
            taxa_abun_data = \
                pd.DataFrame(data=np.vstack((original_taxa_abun_data, residual_function_vs_sample[i, :])),
                             index=np.hstack((original_taxa_abun_data.index.values, "Unknown")),
                             columns=samples)
        else:
            taxa_abun_data = pd.DataFrame(data=original_taxa_abun_data,
                                          index=original_taxa_abun_data.index.values,
                                          columns=samples)

        # now create the full permuted matrices for this function, including
        # the residual (unless removed) and keep the scores of the permuted
        # matrices to subtract later
        permuted_matrices_scores = np.zeros(number_of_permutations)
        permuted_taxa_abundance_matrices = np.zeros((number_of_permutations,
                                                     number_of_samples, num_of_taxa))
        for p in range(number_of_permutations):
            if 'residual_mode' in args.keys() and \
                            args['residual_mode'] != 'remove_residual':
                permuted_taxa_abundance_matrices[p] = \
                    np.hstack((no_residual_permuted_taxa_abundance_matrices[p],
                               np.array([permuted_residual_matrices_function_vs_samples[p, i]]).T))
            else:  # 'remove_residual'
                permuted_taxa_abundance_matrices[p] = \
                    no_residual_permuted_taxa_abundance_matrices[p]

            permuted_taxa_times_weights = np.dot(permuted_taxa_abundance_matrices[p],
                                                 weights_of_this_function)
            permuted_matrices_scores[p] = \
                compute_score(permuted_taxa_times_weights[cases],
                              permuted_taxa_times_weights[controls],
                              args['score_to_compute'], args['max_score_cutoff'])

        # Shapley value analysis with permutations, not of subsets but orderings
        if args['taxa_assessment_method'] == 'multi_taxa':

            print("Computing permuted shapley orderings scores for " +
                  str(number_of_orderings) + " orderings...")

            # permutations are the outer loop and the orderings in the inner loop
            stat_value_for_permutations = np.zeros((number_of_permutations, num_of_taxa))
            for p in range(number_of_permutations):

                absolute_score_value = np.zeros((number_of_orderings, num_of_taxa))
                marginal_score_values = np.zeros((number_of_orderings, num_of_taxa))
                marginal_orderings = np.zeros((number_of_orderings, num_of_taxa))
                # compute shapley value for all orderings
                for o in range(number_of_orderings):
                    current_ordering = random_taxa_orderings[o, :]
                    marginal_orderings[o, :] = current_ordering

                    # now walk across the ordering and compute the score for
                    # every subset in the order
                    for s in range(num_of_taxa):
                        indexing_array = current_ordering[0:(s+1)]
                        non_indexed_array = current_ordering[(s+1):num_of_taxa]

                        curr_matrix_with_permuted_taxa = np.copy(permuted_taxa_abundance_matrices[p])
                        curr_matrix_with_permuted_taxa[:, indexing_array] = \
                            taxa_abun_data.values[indexing_array, :].T

                        if args['normalization_mode'] == 'scale_non_permuted':

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                # residual index is num_of_taxa-1
                                non_residual_indices = (non_indexed_array != (num_of_taxa-1))
                                sum_of_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array[non_residual_indices]], axis=1)
                            else:  # remove_residual or as_taxa
                                sum_of_permuted_taxa_per_sample = np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array], axis=1)

                            target_sum_for_non_permuted = 1.0 - sum_of_permuted_taxa_per_sample

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                # residual index is num_of_taxa-1
                                non_residual_indices = (indexing_array != (num_of_taxa-1))
                                sum_of_non_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, indexing_array[non_residual_indices]], axis=1)
                            else:  # remove_residual or as_taxa
                                sum_of_non_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, indexing_array], axis=1)

                            # if the non-permuted sum to zero, then we must
                            # replace them with a non-zero number, otherwise
                            # they will not be able to be scaled fill the
                            # target sum
                            if sum(sum_of_non_permuted_taxa_per_sample == 0) > 0:

                                for zero_row in np.where(sum_of_non_permuted_taxa_per_sample == 0)[0]:
                                    if 'residual_mode' in args.keys() and \
                                                    args['residual_mode'] == 'as_baseline':
                                        curr_matrix_with_permuted_taxa[zero_row, indexing_array[non_residual_indices]] = 1.0
                                    else:  # remove_residual or as_taxa
                                        curr_matrix_with_permuted_taxa[zero_row, indexing_array] = 1.0

                                if 'residual_mode' in args.keys() and \
                                                args['residual_mode'] == 'as_baseline':
                                    sum_of_non_permuted_taxa_per_sample = \
                                        np.sum(curr_matrix_with_permuted_taxa[:, indexing_array[non_residual_indices]], axis=1)
                                else:  # remove_residual or as_taxa
                                    sum_of_non_permuted_taxa_per_sample = \
                                        np.sum(curr_matrix_with_permuted_taxa[:, indexing_array], axis=1)

                            scaling_factor_for_non_permuted = np.array((1.0 / target_sum_for_non_permuted) *
                                                                       sum_of_non_permuted_taxa_per_sample)

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                curr_matrix_with_permuted_taxa[:, indexing_array[non_residual_indices]] = \
                                    curr_matrix_with_permuted_taxa[:, indexing_array[non_residual_indices]] / \
                                    (np.tile(scaling_factor_for_non_permuted, (len(indexing_array[non_residual_indices]), 1)).T)
                            else:  # remove_residual or as_taxa
                                curr_matrix_with_permuted_taxa[:, indexing_array] = \
                                    curr_matrix_with_permuted_taxa[:, indexing_array] / \
                                    (np.tile(scaling_factor_for_non_permuted, (s+1, 1)).T)


                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                final_sum_of_curr_matrix = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, 0:(num_of_taxa-1)], axis=1)
                            else:  # remove_residual or as_taxa
                                final_sum_of_curr_matrix = \
                                    np.sum(curr_matrix_with_permuted_taxa, axis=1)

                            if max(final_sum_of_curr_matrix) > 1.001 or min(final_sum_of_curr_matrix) < 0.999:
                                print("max: " + str(max(final_sum_of_curr_matrix)) +
                                      " min: " + str(min(final_sum_of_curr_matrix)))
                                sys.exit('Error: sample normalization could not sum up to 1')

                        if args['normalization_mode'] == 'scale_permuted':

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                # residual index is num_of_taxa-1
                                non_residual_indices = (indexing_array != (num_of_taxa-1))
                                sum_of_non_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, indexing_array[non_residual_indices]], axis=1)
                            else:  # remove_residual or as_taxa
                                sum_of_non_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, indexing_array], axis=1)

                            target_sum_for_permuted = 1.0 - sum_of_non_permuted_taxa_per_sample

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                # residual index is num_of_taxa-1
                                non_residual_indices = (non_indexed_array != (num_of_taxa-1))
                                sum_of_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array[non_residual_indices]], axis=1)
                            else:  # remove_residual or as_taxa
                                sum_of_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array], axis=1)

                            # if the permuted sum to zero, then we must replace
                            # them with a non-zero number, otherwise they
                            # will not be able to be scaled fill the target sum
                            if sum(sum_of_permuted_taxa_per_sample == 0) > 0:

                                for zero_row in np.where(sum_of_permuted_taxa_per_sample == 0)[0]:
                                    if 'residual_mode' in args.keys() and args['residual_mode'] == 'as_baseline':
                                        curr_matrix_with_permuted_taxa[zero_row, non_indexed_array[non_residual_indices]] = 1.0
                                    else:  # remove_residual or as_taxa
                                        curr_matrix_with_permuted_taxa[zero_row, non_indexed_array] = 1.0

                                if 'residual_mode' in args.keys() and args['residual_mode'] == 'as_baseline':
                                    sum_of_permuted_taxa_per_sample = np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array[non_residual_indices]], axis=1)
                                else:  # remove_residual or as_taxa
                                    sum_of_permuted_taxa_per_sample = np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array], axis=1)

                            # create scaling factor for permuted. In case the
                            # non-permuted sum to 1 (this could happen do to
                            # numeric reasons since we are normalizing the
                            # taxa by dividing by the sum), so ignore the
                            # warning on division by zero
                            with warnings.catch_warnings():
                                warnings.filterwarnings("ignore", category=RuntimeWarning)
                                scaling_factor_for_permuted = np.array((1.0 / target_sum_for_permuted) *
                                                                       sum_of_permuted_taxa_per_sample)

                            if 'residual_mode' in args.keys() and args['residual_mode'] == 'as_baseline':
                                curr_matrix_with_permuted_taxa[:, non_indexed_array[non_residual_indices]] = \
                                    curr_matrix_with_permuted_taxa[:, non_indexed_array[non_residual_indices]] / \
                                    (np.tile(scaling_factor_for_permuted, (len(non_indexed_array[non_residual_indices]), 1)).T)
                            else:  # remove_residual or as_taxa
                                curr_matrix_with_permuted_taxa[:, non_indexed_array] = \
                                    curr_matrix_with_permuted_taxa[:, non_indexed_array] / \
                                    (np.tile(scaling_factor_for_permuted, (len(non_indexed_array), 1)).T)

                            if 'residual_mode' in args.keys() and args['residual_mode'] == 'as_baseline':
                                final_sum_of_curr_matrix = np.sum(curr_matrix_with_permuted_taxa[:, 0:(num_of_taxa-1)], axis=1)
                            else:  # remove_residual or as_taxa
                                final_sum_of_curr_matrix = np.sum(curr_matrix_with_permuted_taxa, axis=1)


                            if max(final_sum_of_curr_matrix) > 1.001 or min(final_sum_of_curr_matrix) < 0.999:
                                print("max: " + str(max(final_sum_of_curr_matrix)) + " min: " + str(min(final_sum_of_curr_matrix)))
                                sys.exit('Error: sample normalization could not sum up to 1')

                        constructed_abundance_values = np.dot(curr_matrix_with_permuted_taxa,
                                                              weights_of_this_function)
                        absolute_score_value[o, s] = \
                            compute_score(constructed_abundance_values[cases],
                                          constructed_abundance_values[controls],
                                          args['score_to_compute'],
                                          args['max_score_cutoff']) - permuted_matrices_scores[p]
                        if s == 0:  # first value in ordering so no previous subset
                            marginal_score_values[o, s] = absolute_score_value[o, s]
                        else:  # subtract previous score to have only marginal score
                            marginal_score_values[o, s] = absolute_score_value[o, s] - \
                                                          absolute_score_value[o, (s-1)]

                # now compute the average of all marginals per taxa
                # (of all orderings)
                per_taxa_marginals = np.zeros((num_of_taxa, number_of_orderings))
                for t in range(num_of_taxa):
                    curr_taxa_ind = marginal_orderings == t
                    per_taxa_marginals[t, :] = marginal_score_values[curr_taxa_ind]
                    stat_value_for_permutations[p, t] = np.mean(marginal_score_values[curr_taxa_ind])

            print("Done.")


            mean_stat_value[0:num_of_taxa, i] = np.mean(stat_value_for_permutations, axis=0)
            median_stat_value[0:num_of_taxa, i] = np.median(stat_value_for_permutations, axis=0)
            std_stat_value[0:num_of_taxa, i] = np.std(stat_value_for_permutations, axis=0)
            contribution_matrix[0:num_of_taxa, i] = np.mean(stat_value_for_permutations, axis=0)

        else:  # non-shapley permutations based methods

            for j in range(num_of_taxa):  # num_of_taxa

                stat_value_for_permutations = np.zeros(number_of_permutations)

                # Now permute all other taxa
                if args['taxa_assessment_method'] == 'single_taxa':

                    # for every permuted matrix (based on the number of
                    # permutations param given)
                    for p in range(number_of_permutations):

                        # first, copy to a current matrix the permuted one,
                        # and then change only taxon j to be a non-permuted
                        # column
                        curr_matrix_with_permuted_taxa = np.copy(permuted_taxa_abundance_matrices[p])
                        curr_matrix_with_permuted_taxa[:, j] = taxa_abun_data.values[j, :]

                        if args['normalization_mode'] == 'scale_permuted':
                            indexing_array = np.array([j])
                            non_indexed_array = np.delete(range(num_of_taxa), j)

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                # residual index is num_of_taxa-1
                                non_residual_indices = (indexing_array != (num_of_taxa-1))
                                sum_of_non_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, indexing_array[non_residual_indices]], axis=1)
                            else:  # remove_residual or as_taxa
                                sum_of_non_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, indexing_array], axis=1)

                            target_sum_for_permuted = 1.0 - sum_of_non_permuted_taxa_per_sample

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                # residual index is num_of_taxa-1
                                non_residual_indices = (non_indexed_array != (num_of_taxa-1))
                                sum_of_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array[non_residual_indices]], axis=1)
                            else:  # remove_residual or as_taxa
                                sum_of_permuted_taxa_per_sample = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array], axis=1)

                            # if the permuted sum to zero, then we must
                            # replace them with a non-zero number, otherwise
                            # they will not be able to be scaled fill the
                            # target sum
                            if sum(sum_of_permuted_taxa_per_sample == 0) > 0:

                                for zero_row in np.where(sum_of_permuted_taxa_per_sample == 0)[0]:
                                    if 'residual_mode' in args.keys() and \
                                                    args['residual_mode'] == 'as_baseline':
                                        curr_matrix_with_permuted_taxa[zero_row, non_indexed_array[non_residual_indices]] = 1.0
                                    else:  # remove_residual or as_taxa
                                        curr_matrix_with_permuted_taxa[zero_row, non_indexed_array] = 1.0

                                if 'residual_mode' in args.keys() and \
                                                args['residual_mode'] == 'as_baseline':
                                    sum_of_permuted_taxa_per_sample = \
                                        np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array[non_residual_indices]], axis=1)
                                else:  # remove_residual or as_taxa
                                    sum_of_permuted_taxa_per_sample = \
                                        np.sum(curr_matrix_with_permuted_taxa[:, non_indexed_array], axis=1)


                            scaling_factor_for_permuted = np.array((1.0 / target_sum_for_permuted) *
                                                                   sum_of_permuted_taxa_per_sample)

                            curr_matrix_with_permuted_taxa[:, non_indexed_array] = \
                                curr_matrix_with_permuted_taxa[:, non_indexed_array] / \
                                (np.tile(scaling_factor_for_permuted, (len(non_indexed_array), 1)).T)

                            if 'residual_mode' in args.keys() and \
                                            args['residual_mode'] == 'as_baseline':
                                final_sum_of_curr_matrix = \
                                    np.sum(curr_matrix_with_permuted_taxa[:, 0:(num_of_taxa-1)], axis=1)
                            else:  # remove_residual or as_taxa
                                final_sum_of_curr_matrix = \
                                    np.sum(curr_matrix_with_permuted_taxa, axis=1)

                            if max(final_sum_of_curr_matrix) > 1.001 or min(final_sum_of_curr_matrix) < 0.999:
                                print("max: " + str(max(final_sum_of_curr_matrix)) + " min: " +
                                      str(min(final_sum_of_curr_matrix)))
                                sys.exit('Error: sample normalization could not sum up to 1')

                        constructed_abundance_values = \
                            np.dot(curr_matrix_with_permuted_taxa, weights_of_this_function)
                        stat_value_for_permutations[p] = \
                            compute_score(constructed_abundance_values[cases],
                                          constructed_abundance_values[controls],
                                          args['score_to_compute'],
                                          args['max_score_cutoff']) - permuted_matrices_scores[p]

                    mean_stat_value[j, i] = np.mean(stat_value_for_permutations)
                    median_stat_value[j, i] = np.median(stat_value_for_permutations)
                    std_stat_value[j, i] = np.std(stat_value_for_permutations)
                    contribution_matrix[j, i] = np.mean(stat_value_for_permutations)

                # For testing purposes, permute all taxa
                elif args['taxa_assessment_method'] == 'permute_all':

                    for p in range(number_of_permutations):
                        curr_matrix_with_permuted_taxa = \
                            np.copy(permuted_taxa_abundance_matrices[p])
                        constructed_abundance_values = \
                            np.dot(curr_matrix_with_permuted_taxa, weights_of_this_function)
                        stat_value_for_permutations[p] = \
                            compute_score(constructed_abundance_values[cases],
                                          constructed_abundance_values[controls],
                                          args['score_to_compute'],
                                          args['max_score_cutoff']) - permuted_matrices_scores[p]

                    mean_stat_value[j, i] = np.mean(stat_value_for_permutations)
                    median_stat_value[j, i] = np.median(stat_value_for_permutations)
                    std_stat_value[j, i] = np.std(stat_value_for_permutations)
                    contribution_matrix[j, i] = np.mean(stat_value_for_permutations)

        end = time.time()
        print(str(i) + ":" + function_abun_data.index.values[i] + " took " +
              str(end - start) + " seconds to run.")
        if 'write_log' in args.keys() and args['write_log']:
            with open(args['output_pref'] + '_STAT_run_log' + output_suffix, 'a') as f:
                f.write(str(i) + ":" + function_abun_data.index.values[i] +
                        " took " + str(end - start) + " seconds to run." + "\n")

    ###########################################################################
    # WRITE OUTPUT
    ###########################################################################

    print("Writing output...")

    # write the FUNCTION differential abundance results into file
    with open(args['output_pref'] + '_STAT_DA_function' + output_suffix, 'w') as f:
        f.write("# " + sys.argv[0] + " " + str(args) + '\n')
    functions_da_scores.index.name = 'Function'
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        functions_da_scores.to_csv(args['output_pref'] + '_STAT_DA_function' +
                                   output_suffix, sep='\t', na_rep=args['na_rep'],
                                   mode='a')

    # write the TAXA differential abundance results into file
    with open(args['output_pref'] + '_STAT_DA_taxa' + output_suffix, 'w') as f:
        f.write("# " + sys.argv[0] + " " + str(args) + '\n')
    taxa_diff_abun_scores.index.name = 'Taxa'
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        taxa_diff_abun_scores.to_csv(args['output_pref'] + '_STAT_DA_taxa' +
                                     output_suffix, sep='\t', na_rep=args['na_rep'],
                                     mode='a')

    # write the expected functional abundance into file
    cont_pd = pd.DataFrame(data=predicted_function_abundance.T,
                           index=function_abun_data.index.values,
                           columns=function_abun_data.columns.values)
    cont_pd.index.name = 'KO'
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        cont_pd.to_csv(args['output_pref'] +
                       '_STAT_predicted_function_abundance' + output_suffix,
                       sep='\t', na_rep=args['na_rep'])

    # write the residual of functional abundance into file
    cont_pd = pd.DataFrame(data=residual_function_vs_sample,
                           index=function_abun_data.index.values,
                           columns=function_abun_data.columns.values)
    cont_pd.index.name = 'KO'
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        cont_pd.to_csv(args['output_pref'] +
                       '_STAT_residual_function_abundance' + output_suffix,
                       sep='\t', na_rep=args['na_rep'])

    if num_of_da_functions > 0:
        # print out the main output file, detailing the contribution of each
        # taxon to each functional shift including the distinction between 4
        #  modes of contribution
        main_output_matrix = contribution_matrix.astype(object)
        for i in range(num_of_da_functions):
            for t in range(num_of_taxa):
                if taxa_diff_abun_scores.values[t, 2] > 0:
                    if contribution_matrix[t, i] > 0:
                        main_output_matrix[t, i] = "a:" + str(main_output_matrix[t, i])
                    else:
                        main_output_matrix[t, i] = "b:" + str(main_output_matrix[t, i])
                else:
                    if contribution_matrix[t, i] > 0:
                        main_output_matrix[t, i] = "c:" + str(main_output_matrix[t, i])
                    else:
                        main_output_matrix[t, i] = "d:" + str(main_output_matrix[t, i])

        cont_pd = pd.DataFrame(data=main_output_matrix,
                               index=taxa_abun_data.index.values,
                               columns=function_abun_data.index.values)
        cont_pd.index.name = 'Taxa'
        with open(args['output_pref'] + '_main_output' + output_suffix, 'w') as f:
            f.write("# RUN PARAMETERS: " + '\n')
            f.write("# --------------- " + '\n')
            f.write("# " + sys.argv[0] + " " + str(args) + '\n')
            f.write("#" + '\n')
            f.write("# LEGEND: " + '\n')
            f.write("# ------- " + '\n')
            f.write("# a:case-associated and driving case-enrichment; " + '\n')
            f.write("# b:case-associated and attenuating case-enrichment; " + '\n')
            f.write("# c:control-associated and driving case-enrichment; " + '\n')
            f.write("# d:control-associated and attenuating case-enrichment; " + '\n')
            f.write("#" + '\n')
            f.write("# TAXON-LEVEL CONTRIBUTION VALUES: " + '\n')
            f.write("# -------------------------------- " + '\n')
            f.write("#" + '\n')
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_main_output' + output_suffix,
                           sep='\t', na_rep=args['na_rep'], encoding='utf-8', mode='a')

        # print out the supporting contribution value of each taxon to each functional shift
        cont_pd = pd.DataFrame(data=contribution_matrix, index=taxa_abun_data.index.values,
                               columns=function_abun_data.index.values)
        cont_pd.index.name = 'Taxa'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_STAT_taxa_contributions' +
                           output_suffix, sep='\t', na_rep=args['na_rep'])

        cont_pd = pd.DataFrame(data=original_stat_value, index=function_abun_data.index.values,
                               columns=[args['score_to_compute']])
        cont_pd.index.name = 'KO'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_STAT_original_value' +
                           output_suffix, sep='\t', na_rep=args['na_rep'])

        cont_pd = pd.DataFrame(data=predicted_da_stat_value, index=function_abun_data.index.values,
                               columns=[args['score_to_compute']])
        cont_pd.index.name = 'KO'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_STAT_predicted_DA_value' +
                           output_suffix, sep='\t', na_rep=args['na_rep'])

        cont_pd = pd.DataFrame(data=predicted_function_agreement, index=function_abun_data.index.values,
                               columns=np.array(("R^2", "PearsonCorr", "PearsonPval", "SpearmanCorr",
                                                 "SpearmanPval", "MeanAbsDiff", "StdAbsDiff")))
        cont_pd.index.name = 'KO'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_STAT_predicted_function_agreement' +
                           output_suffix, sep='\t', na_rep=args['na_rep'])

        cont_pd = pd.DataFrame(data=mean_stat_value, index=taxa_abun_data.index.values,
                               columns=function_abun_data.index.values)
        cont_pd.index.name = 'Taxa'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_STAT_mean_stat' +
                           output_suffix, sep='\t', na_rep=args['na_rep'])

        cont_pd = pd.DataFrame(data=median_stat_value, index=taxa_abun_data.index.values,
                               columns=function_abun_data.index.values)
        cont_pd.index.name = 'Taxa'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_STAT_median_stat' +
                           output_suffix, sep='\t', na_rep=args['na_rep'])

        cont_pd = pd.DataFrame(data=std_stat_value, index=taxa_abun_data.index.values,
                               columns=function_abun_data.index.values)
        cont_pd.index.name = 'Taxa'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            cont_pd.to_csv(args['output_pref'] + '_STAT_std_stat' +
                           output_suffix, sep='\t', na_rep=args['na_rep'])

        # write the inferred (or original if no inference was done) genome
        # content in each taxa
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            original_taxa_to_function_data.to_csv(args['output_pref'] +
                                                  '_STAT_taxa_learned_copy_num' +
                                                  output_suffix, sep='\t', na_rep='None')

        # print out R^2 and learned taxa-to-function values, from either real
        # inference or given data
        if 'apply_inference' in args.keys() and args['apply_inference']:
            # write the test rsqr and other stats results into file
            global_cv_pd = pd.DataFrame(data=all_functions_global_cv_test_stats,
                                        index=functions_to_infer,
                                        columns=['Global_Test_RSQR', 'Global_Test_Pearson',
                                                 'Global_Test_Spearman'])
            global_cv_pd.index.name = 'KO'
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=DeprecationWarning)
                global_cv_pd.to_csv(args['output_pref'] + '_STAT_taxa_learning_rsqr' +
                                    output_suffix, sep='\t', na_rep='None')

        else:
            # copy the function agreement (only the R^2 part) as the R^2 of
            # these functions
            global_cv_pd = pd.DataFrame(data=predicted_function_agreement[:, [0, 1, 3]],
                                        index=functions,
                                        columns=['Global_Test_RSQR', 'Global_Test_Pearson',
                                                 'Global_Test_Spearman'])
            global_cv_pd.index.name = 'KO'
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=DeprecationWarning)
                global_cv_pd.to_csv(args['output_pref'] + '_STAT_taxa_learning_rsqr' +
                                    output_suffix, sep='\t', na_rep='None')

        print("Done.")

    else:
        print("Note that since there are no differentially abundant functions, "
              "there is no real output...")

    if 'write_log' in args.keys() and args['write_log']:
        with open(args['output_pref'] + '_STAT_run_log' + output_suffix, 'a') as f:
            f.write("----------------------------------------------" + "\n")
            f.write("Program completed successfully with no errors." + "\n")
            f.write("----------------------------------------------" + "\n")

    return

###############################################################################

if __name__ == "__main__":
    # get options from user
    parser = \
        argparse.ArgumentParser(description='Quantify the individual '
                                            'contributions of taxa to shifts '
                                            'observed in functional '
                                            'composition across different '
                                            'sample sets')

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
                        help='Input file of label assignment for the two sample '
                             'sets being compared', default=None)

    # Optional arguments:

    parser.add_argument('-gc', '--genomic_content_of_taxa',
                        dest='taxa_to_function_file',
                        help='Input file of genomic content of each taxa',
                        default=None)

    parser.add_argument('-inf', '--perform_inference_of_genomic_content',
                        dest='apply_inference',
                        help='Defines if genome content is inferred '
                             '(either de-novo or prior-based if genomic '
                             'content is also given)', action='store_true')

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
                        help='pathways or modules mapping file (default: '
                             'use internal KEGG file)', default=None)

    parser.add_argument('-perform_inference_on_ko_level',
                        dest='perform_inference_on_ko_level',
                        help='Indicates to perform the inference on the KO '
                             'level (default: use the mapped functional level, '
                             'e.g., pathway)', action='store_true')

    parser.add_argument('-mult_hyp', '--multiple_hypothesis_correction',
                        dest='multiple_hypothesis_correction',
                        help='Multiple hypothesis correction for functional '
                             'enrichment (default: FDR-0.05)',
                        default='FDR-0.05',
                        choices=['Bonf', 'FDR-0.01', 'FDR-0.05',
                                 'FDR-0.1', 'none'])

    parser.add_argument('-max_func', '--maximum_functions_to_analyze',
                        dest='max_da_functions_cases_controls',
                        help='Maximum number of enriched functions to '
                             'consider (default: All)', default=None)

    parser.add_argument('-assessment', '--taxa_assessment_method',
                        dest='taxa_assessment_method',
                        help='The method used when assessing taxa to '
                             'compute individual contributions (default: '
                             'multi_taxa)',
                        default='multi_taxa', choices=['single_taxa',
                                                       'multi_taxa'])

    parser.add_argument('-score', '--score_to_compute',
                        dest='score_to_compute',
                        help='The enrichment score to compute for each '
                             'function (default: wilcoxon)',
                        default='wilcoxon', choices=['t_test', 'mean_diff',
                                                     'median_diff', 'wilcoxon',
                                                     'log_mean_ratio'])

    parser.add_argument('-max_score', '--max_score_cutoff',
                        dest='max_score_cutoff',
                        help='The maximum score cutoff (for example, '
                             'when dividing by zero) (default: 100)',
                        default='100')

    parser.add_argument('-na_rep', dest='na_rep',
                        help='How to represent NAs in the output (default: NA)',
                        default='NA')

    parser.add_argument('-number_of_permutations',
                        dest='number_of_permutations',
                        help='number of permutations (default: 100)',
                        default='100')

    parser.add_argument('-number_of_shapley_orderings_per_taxa',
                        dest='number_of_shapley_orderings_per_taxa',
                        help='number of shapley orderings per taxa '
                             '(default: 5)', default='5')

    # DEPRECATED:
    # parser.add_argument('-use_gc_as_prior',
    # '--use_genomic_content_of_taxa_as_prior',
    # dest='use_t2f_as_prior',
    #                    help='Learn the taxa copy number of each function,
    # using the given genomic content data as prior (default: False)',
    # action='store_true')
    # parser.add_argument('-residual_mode', dest='residual_mode',
    # choices=['as_taxa', 'remove_residual', 'as_baseline'],
    #                    help='How to treat the residual of the functional
    # abundance profile (default: remove_residual)', default='remove_residual')
    # parser.add_argument('-normalization_mode', dest='normalization_mode',
    # choices=['none', 'scale_non_permuted', 'scale_permuted'],
    #                    help='How to normalize the sample after permuting
    # taxa (default: scale_permuted)', default='scale_permuted')
    # parser.add_argument('-permutation_mode', dest='permutation_mode',
    # choices=['independent', 'blocks'],
    #                    help='How to permute the taxa across samples
    # (default: blocks)', default='blocks')

    parser.add_argument('-en', '--enrichment_results', dest='da_result_file',
                        help='Pre-computed functional enrichment results from '
                             'the compute_differential_abundance.py script '
                             '(default: None)', default=None)

    parser.add_argument('-single_function_filter',
                        dest='single_function_filter',
                        help='Limit analysis only to this single '
                             'function (default: None)', default=None)

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





