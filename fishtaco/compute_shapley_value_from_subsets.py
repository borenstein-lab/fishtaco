"""
This function computes the Shapley value for each taxon from the given
shift scores given for many subsets of taxa

Parameters
    ----------

    args: dictionary
        a dictionary containing the function parameters
        args.['shapley_input_prefix']: The prefix of the input file of shapley
        values
        args.['shapley_input_matrix_name']: The input file matrix name of
        shapley values
        args.['shapley_input_functions']: The different functions to compute
        for
        args.['subset_assessment_method']: The method used when assessing
        subsets to compute Shapley scores
        args.['subset_score_method']: The score used to compute for Shapley
        values each subset
        args.['shapley_input_suffix']: The suffix of the input file of shapley
        values
        args.['marginals_output_prefix']: The prefix of the output file for
        marginal shapley values
        args.['shapley_output_matrix_name']: The output file matrix name of
        shapley values
        args.['marginals_output_suffix']: The suffix of the output file for
        marginal shapley values
        args.['taxa_DA_input_file_suffix']: The input file suffix of the
        differential abundance of taxa
        args.['output_per_function']: Output Shapley marginal values per
        function

"""

import numpy as np
import pandas as pd
import argparse
import re
import sys
import os

__author__ = 'Ohad Manor'
__email__ = 'omanor@gmail.com'
__status__ = "Development"


###############################################################################
# COMMENT COUNTER FUNCTION - counts the number of comments in the file
###############################################################################
def n_comments(fn, comment):
    with open(fn, 'r') as f:
        n_lines = 0
        pattern = re.compile("^\s*{0}".format(comment))
        for l in f:
            if pattern.search(l) is None:
                break
            else:
                n_lines += 1
    return n_lines


###############################################################################
# MAIN FUNCTION
###############################################################################
def main(args):

    print("Given parameters: ", args)
    pd.set_option('display.max_rows', 500)

    given_functions = args['shapley_input_functions'].split(",")
    num_of_functions = len(given_functions)
    print("list of " + str(num_of_functions) + " given functions: " +
          str(given_functions))
    output_file_all_functions = args['marginals_output_prefix'] + "_" +\
                                args['subset_score_method'] + "_" +\
                                args['shapley_output_matrix_name'] + "_" +\
                                args['subset_assessment_method'] +\
                                args['shapley_input_suffix']


    ########################################################
    # READ THE TAXA DA FILE, NOTE THAT WE ARE SKIPPING THE
    # FIRST LINES THAT BEGIN WITH A COMMENT CHAR "#"
    ########################################################

    if args['taxa_DA_input_file_suffix'] is not None:
        taxa_da_file = args['shapley_input_prefix'] + "_" + \
                       args['subset_score_method'] + "_" + \
                       args['subset_assessment_method'] + "_" + \
                       args['taxa_DA_input_file_suffix']
        if not os.path.isfile(taxa_da_file):
            sys.exit('Error: Input file "' + taxa_da_file + '" does not exist')
        taxa_da_pd = pd.read_table(taxa_da_file,
                                   skiprows=n_comments(taxa_da_file, '#'),
                                   index_col=0)
        num_of_taxa = taxa_da_pd.index.values.shape[0]
    else:
        sys.exit('Error: No taxa DA input given to script')

    ########################################################
    # COMPUTE THE SHAPLEY VALUES FOR EACH GIVEN FUNCTION
    ########################################################

    median_of_marginals = np.zeros((num_of_taxa, num_of_functions))

    for f in range(num_of_functions):

        print("current function: " + given_functions[f])

        ########################################################
        # READ THE SUBSET SCORE VALUES FILE, NOTE THAT WE ARE SKIPPING THE
        # FIRST LINES THAT BEGIN WITH A COMMENT CHAR "#"
        ########################################################

        if args['shapley_input_prefix'] is not None:
            input_file = args['shapley_input_prefix'] + "_" + \
                         args['subset_score_method'] + "_" + \
                         args['shapley_input_matrix_name'] + "_" + \
                         args['subset_assessment_method'] + "_" + \
                         given_functions[f] + args['shapley_input_suffix']
            if not os.path.isfile(input_file):
                sys.exit('Error: Input file "' + input_file +
                         '" does not exist')
            subset_data = pd.read_table(input_file,
                                        skiprows=n_comments(input_file, '#'),
                                        index_col=0)
        else:
            sys.exit('Error: No shapley scores input given to script')

        if args['marginals_output_prefix'] is not None:
            output_file_per_function = \
                args['marginals_output_prefix'] + "_" + \
                args['subset_score_method'] + "_" + \
                args['shapley_output_matrix_name'] + "_" + \
                args['subset_assessment_method'] + "_" + \
                given_functions[f] + args['shapley_input_suffix']
        else:
            sys.exit('Error: No shapley output prefix given to script')

        ########################################################
        # ANALYZE THE FILE
        ########################################################

        if args['subset_assessment_method'] == 'shapley' or \
                        args['subset_assessment_method'] == 'permuted_shapley':
            num_of_taxa = subset_data.shape[1] - 1

            size_of_subset = np.sum(subset_data.values[:, 0:num_of_taxa],
                                    axis=1)

            subset_sizes = [1, 2, 3, 4, 5, num_of_taxa-3, num_of_taxa-2,
                            num_of_taxa-1, num_of_taxa]

            average_of_marginals = np.zeros((num_of_taxa, len(subset_sizes)))

            # i loops over all the taxa range(num_of_taxa)
            for i in range(num_of_taxa):

                print("i: " + str(i))

                subset_counter = 0

                index_without_curr = np.delete(np.arange(num_of_taxa), i)

                for s in subset_sizes:

                    if s == 1:
                        locations_of_size_s = (subset_data.values[:, i] == 1) \
                                              & (size_of_subset == 1)
                        average_of_marginals[i, subset_counter] = \
                            subset_data.values[locations_of_size_s,
                                               num_of_taxa]

                    elif s == num_of_taxa:
                        score_for_size_s = subset_data.values[size_of_subset ==
                                                              num_of_taxa, num_of_taxa][0]
                        locations_of_size_s_minus_1 = (subset_data.values[:, i] == 0) \
                                                      & (size_of_subset == s-1)
                        score_for_size_s_minus_1 = \
                            subset_data.values[locations_of_size_s_minus_1,
                                               num_of_taxa][0]
                        average_of_marginals[i, subset_counter] = \
                            score_for_size_s - score_for_size_s_minus_1

                    elif s == 2:
                        locations_of_size_s = (subset_data.values[:, i] == 1) \
                                              & (size_of_subset == s)
                        subset_nonzero_index_of_size_s = \
                            np.reshape(np.nonzero(subset_data.values[np.ix_(locations_of_size_s,
                                                                            index_without_curr)])[1],
                                       (np.sum(locations_of_size_s), s-1))
                        string_rep_of_nonzero_of_size_s = \
                            np.char.mod('%d', subset_nonzero_index_of_size_s)[:, 0]
                        scores_for_subsets_of_size_s = \
                            subset_data.values[locations_of_size_s, num_of_taxa]
                        df_of_size_s = \
                            pd.DataFrame(data=scores_for_subsets_of_size_s,
                                         index=string_rep_of_nonzero_of_size_s,
                                         columns=['size_s'])

                        # only for s == 2
                        locations_of_size_s_minus_1 = (subset_data.values[:, i] == 0) \
                                                      & (size_of_subset == s-1)
                        subset_nonzero_index_of_size_s_minus_1 = \
                            np.reshape(np.nonzero(subset_data.values[np.ix_(locations_of_size_s_minus_1,
                                                                            index_without_curr)])[1],
                                       (np.sum(locations_of_size_s_minus_1), s-1))
                        scores_for_subsets_of_size_s_minus_1 = \
                            subset_data.values[locations_of_size_s_minus_1,
                                               num_of_taxa]
                        string_rep_of_nonzero_of_size_s_minus_1 = \
                            np.char.mod('%d', subset_nonzero_index_of_size_s_minus_1)[:, 0]
                        df_of_size_s_minus_1 = \
                            pd.DataFrame(data=scores_for_subsets_of_size_s_minus_1,
                                         index=string_rep_of_nonzero_of_size_s_minus_1,
                                         columns=['size_s_minus_1'])

                        average_of_marginals[i, subset_counter] = \
                            np.median(df_of_size_s.join(df_of_size_s_minus_1)['size_s'] -
                                      df_of_size_s.join(df_of_size_s_minus_1)['size_s_minus_1'])

                    else:
                        # first for size s

                        locations_of_size_s = (subset_data.values[:, i] == 1) & \
                                              (size_of_subset == s)
                        delim_array = np.empty(np.sum(locations_of_size_s),
                                               dtype='<U3')
                        delim_array[:] = '_'

                        if s < (num_of_taxa / 2):  # small number of 1's
                            subset_nonzero_index_of_size_s = \
                                np.reshape(np.nonzero(subset_data.values[np.ix_(locations_of_size_s,
                                                                                index_without_curr)])[1],
                                           (np.sum(locations_of_size_s), s-1))
                            string_rep_of_nonzero_of_size_s = \
                                np.char.mod('%d', subset_nonzero_index_of_size_s[:, 0])
                            for j in np.arange(1, s-1):
                                string_rep_of_nonzero_of_size_s = \
                                    np.core.defchararray.add(np.core.defchararray.add(string_rep_of_nonzero_of_size_s,
                                                                                      delim_array),
                                                             np.char.mod('%d', subset_nonzero_index_of_size_s[:, j]))
                        else:  # small number of 0's
                            subset_nonzero_index_of_size_s = \
                                np.reshape(np.nonzero(1 - subset_data.values[np.ix_(locations_of_size_s,
                                                                                    index_without_curr)])[1],
                                           (np.sum(locations_of_size_s), num_of_taxa - s))
                            string_rep_of_nonzero_of_size_s = \
                                np.char.mod('%d', subset_nonzero_index_of_size_s[:, 0])
                            for j in np.arange(1, num_of_taxa - s):
                                string_rep_of_nonzero_of_size_s = \
                                    np.core.defchararray.add(np.core.defchararray.add(string_rep_of_nonzero_of_size_s, delim_array),
                                                             np.char.mod('%d', subset_nonzero_index_of_size_s[:, j]))

                        scores_for_subsets_of_size_s = \
                            subset_data.values[locations_of_size_s, num_of_taxa]
                        df_of_size_s = pd.DataFrame(data=scores_for_subsets_of_size_s,
                                                    index=string_rep_of_nonzero_of_size_s,
                                                    columns=['size_s'])

                        # now for size of subset s-1
                        locations_of_size_s_minus_1 = (size_of_subset == s-1)
                        is_curr_taxa_absent = \
                            np.sum(subset_data.values[np.ix_(locations_of_size_s_minus_1,
                                                             index_without_curr)], axis=1) == s-1
                        locations_of_size_s_minus_1[locations_of_size_s_minus_1] = \
                            is_curr_taxa_absent
                        delim_array = np.empty(np.sum(locations_of_size_s_minus_1),
                                               dtype='<U3')
                        delim_array[:] = '_'

                        if s < (num_of_taxa / 2):  # small number of 1's
                            subset_nonzero_index_of_size_s_minus_1 = \
                                np.reshape(np.nonzero(subset_data.values[np.ix_(locations_of_size_s_minus_1,
                                                                                index_without_curr)])[1],
                                           (np.sum(locations_of_size_s_minus_1), s-1))
                            string_rep_of_nonzero_of_size_s_minus_1 = \
                                np.char.mod('%d', subset_nonzero_index_of_size_s_minus_1[:, 0])
                            for j in np.arange(1, s-1):
                                string_rep_of_nonzero_of_size_s_minus_1 = \
                                    np.core.defchararray.add(np.core.defchararray.add(string_rep_of_nonzero_of_size_s_minus_1, delim_array),
                                                             np.char.mod('%d', subset_nonzero_index_of_size_s_minus_1[:, j]))
                        else:  # small number of 0's
                            subset_nonzero_index_of_size_s_minus_1 = \
                                np.reshape(np.nonzero(1 - subset_data.values[np.ix_(locations_of_size_s_minus_1, index_without_curr)])[1],
                                           (np.sum(locations_of_size_s_minus_1), num_of_taxa - s))
                            string_rep_of_nonzero_of_size_s_minus_1 = \
                                np.char.mod('%d', subset_nonzero_index_of_size_s_minus_1[:, 0])
                            for j in np.arange(1, num_of_taxa - s):
                                string_rep_of_nonzero_of_size_s_minus_1 = \
                                    np.core.defchararray.add(np.core.defchararray.add(string_rep_of_nonzero_of_size_s_minus_1, delim_array),
                                                             np.char.mod('%d', subset_nonzero_index_of_size_s_minus_1[:, j]))


                        scores_for_subsets_of_size_s_minus_1 = \
                            subset_data.values[locations_of_size_s_minus_1, num_of_taxa]
                        df_of_size_s_minus_1 = \
                            pd.DataFrame(data=scores_for_subsets_of_size_s_minus_1,
                                         index=string_rep_of_nonzero_of_size_s_minus_1,
                                         columns=['size_s_minus_1'])

                        average_of_marginals[i, subset_counter] = \
                            np.median(df_of_size_s.join(df_of_size_s_minus_1, how='inner')['size_s'] -
                                      df_of_size_s.join(df_of_size_s_minus_1, how='inner')['size_s_minus_1'])

                    subset_counter += 1

            if output_file_per_function:
                df_marginals = pd.DataFrame(data=average_of_marginals,
                                            index=taxa_da_pd.index.values,
                                            columns=subset_sizes)
                df_marginals.index.name = 'Taxa'
                df_marginals.to_csv(output_file_per_function, sep='\t',
                                    na_rep='NA')

            # keep the median of marginals for all functions
            median_of_marginals[:, f] = np.median(average_of_marginals, axis=1)

        # permuted_shapley_orderings go over all orderings and turn them into
        # marginals, and then compute median for each taxa
        else:
            num_of_taxa = subset_data.shape[1] // 2
            num_of_orderings = subset_data.shape[0]
            print("num of taxa: " + str(num_of_taxa) + ", num of orderings: " +
                  str(num_of_orderings))
            print("Computing marginals from orderings...")
            marginal_score_values = np.zeros((num_of_orderings,num_of_taxa))

            for o in range(num_of_orderings):
                marginal_score_values[o, 0] = subset_data.values[o, num_of_taxa]
                for t in np.arange((num_of_taxa+1), (num_of_taxa*2)):
                   marginal_score_values[o, t-num_of_taxa] = \
                       subset_data.values[o, t] - subset_data.values[o, (t-1)]

            # compute the average of all marginals per taxa (of all orderings)
            per_taxa_marginals = np.zeros((num_of_taxa, num_of_orderings))
            taxa_orderings = subset_data.values[:, 0:num_of_taxa]
            for t in range(num_of_taxa):
                curr_taxa_ind = taxa_orderings == t
                per_taxa_marginals[t, :] = marginal_score_values[curr_taxa_ind]

            if output_file_per_function:
                df_marginals = pd.DataFrame(data=per_taxa_marginals,
                                            index=taxa_da_pd.index.values,
                                            columns=np.arange(1,num_of_orderings+1))
                df_marginals.index.name = 'Taxa'
                df_marginals.to_csv(output_file_per_function, sep='\t',
                                    na_rep='NA')

            # keep the median of marginals for all functions
            median_of_marginals[:, f] = np.median(per_taxa_marginals, axis=1)
            print("Done.")

    print("Writing output...")

    df_median = pd.DataFrame(data=median_of_marginals,
                             index=taxa_da_pd.index.values,
                             columns=given_functions)
    df_median.index.name = 'Taxa'
    df_median.to_csv(output_file_all_functions, sep='\t', na_rep='NA')

    print("Done.")

###############################################################################

if __name__ == "__main__":

    # get options from user
    parser = \
        argparse.ArgumentParser(description='Estimate the shapley value for '
                                            'each taxa from subset scores '
                                            'created by the script '
                                            'compute_contribution_to_DA.py')

    parser.add_argument('-ip', '--shapley_input_prefix',
                        dest='shapley_input_prefix',
                        help='The prefix of the input file of shapley values',
                        default=None)

    parser.add_argument('-im', '--shapley_input_matrix_name',
                        dest='shapley_input_matrix_name',
                        help='The input file matrix name of shapley values ('
                             'default: taxa_contributions_PERMUTATION)',
                        default="taxa_contributions_PERMUTATION")

    parser.add_argument('-if', '--shapley_input_functions',
                        dest='shapley_input_functions',
                        help='The different functions to compute for',
                        default=None)

    parser.add_argument('-assessment', '--subset_assessment_method',
                        dest='subset_assessment_method',
                        help='The method used when assessing subsets to '
                             'compute Shapley scores (default: '
                             'permuted_shapley)',
                        default='permuted_shapley',
                        choices=['shapley', 'permuted_shapley',
                                 'permuted_shapley_orderings'])

    parser.add_argument('-score', '--subset_score_method',
                        dest='subset_score_method',
                        help='The score used to compute for Shapley values '
                             'each subset (default: t-test)',
                        default='t_test',
                        choices=['t_test', 'mean_diff', 'median_diff',
                                 'wilcoxon', 'mean_ratio', 'log_mean_ratio'])

    parser.add_argument('-is', '--shapley_input_suffix',
                        dest='shapley_input_suffix',
                        help='The suffix of the input file of shapley '
                             'values (default: .tab)', default=".tab")

    parser.add_argument('-op', '--marginals_output_prefix',
                        dest='marginals_output_prefix',
                        help='The prefix of the output file for marginal '
                             'shapley values (default: shapley)',
                        default="shapley")

    parser.add_argument('-om', '--shapley_output_matrix_name',
                        dest='shapley_output_matrix_name',
                        help='The output file matrix name of shapley values '
                             '(default: taxa_contributions_MARGINALS_SHAPLEY)',
                        default="taxa_contributions_MARGINALS_SHAPLEY")

    parser.add_argument('-os', '--marginals_output_suffix',
                        dest='marginals_output_suffix',
                        help='The suffix of the output file for marginal '
                             'shapley values (default: .tab)', default=".tab")

    parser.add_argument('-ts', '--taxa_DA_input_file_suffix',
                        dest='taxa_DA_input_file_suffix',
                        help='The input file suffix of the differential '
                             'abundance of taxa (default: DA_taxa.tab)',
                        default="DA_taxa.tab")

    parser.add_argument('-opf', '--output_per_function',
                        dest='output_per_function',
                        help='Output Shapley marginal values per function ('
                             'default: false)', action='store_true')

    given_args = parser.parse_args()
    main(vars(given_args))
