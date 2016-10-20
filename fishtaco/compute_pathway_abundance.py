"""
This function aggregates the abundance of KOs into pathways or modules
abundances and writes the output to the given file

Parameters
----------

    args: dictionary
        a dictionary containing the function parameters
        args.['ko_abun_file']: Input file of ko abundance per sample
        args.['ko_to_pathway_file']: Input file of mappingfrom ko to pathway
        args.['output_file']: Output file for resulting pathway abundance
        args.['output_counts_file']: Output file for number of KOs mapped to
        each pathway
        args.['mapping_method']: Method to map KOs to Pathway
        args.['compute_method']: Method to compute pathway abundance from
        mapped KOs
        args.['transpose_ko_abundance']: Transpose the ko abundance matrix given
        args.['transpose_output']: Transpose the output pathway abundance matrix
        args.['verbose']: Increase verbosity of module

"""

# to comply with both Py2 and Py3
from __future__ import absolute_import, division, print_function

import argparse
import numpy as np
import pandas as pd
import os
import sys
import warnings

__author__ = 'Ohad Manor'
__email__ = 'omanor@gmail.com'
__status__ = "Development"

###############################################################################
# MAIN FUNCTION
###############################################################################
def main(args):

    if 'verbose' in args.keys() and args['verbose']:
        print("Given parameters: ", args)

    ###########################################################################
    # INPUT
    ###########################################################################
    print("Reading files...")

    if 'ko_abun_file' in args.keys() and args['ko_abun_file'] is not None:
        if not os.path.isfile(args['ko_abun_file']):
            sys.exit('Error: Input file "' +
                     args['ko_abun_file'] + '" does not exist')
        ko_abun_data = pd.read_table(args['ko_abun_file'], index_col=0,
                                     dtype={0: str})

    elif 'ko_abun_pd' in args.keys():
        ko_abun_data = args['ko_abun_pd']

    else:
        sys.exit('Error: No input ko abundance file given to script')

    if args['transpose_ko_abundance']:
        ko_abun_data = ko_abun_data.T

    if 'output_file' in args.keys() and args['ko_to_pathway_file'] is not None:
        if not os.path.isfile(args['ko_to_pathway_file']):
            sys.exit('Error: Input file "' +
                     args['ko_to_pathway_file'] + '" does not exist')
        ko_to_pathway_data = pd.read_table(args['ko_to_pathway_file'],
                                           index_col=0, dtype={0: str})

    else:
        if args['ko_to_pathway_pd'] is not None:
            ko_to_pathway_data = args['ko_to_pathway_pd']
        else:
            sys.exit('Error: No input ko to pathway file given to script')

    print("Done.")

    ###########################################################################
    # FILTER OUT KOs THAT ARE NOT PART OF ANY PATHWAY
    ###########################################################################
    ko = np.sort(np.intersect1d(ko_abun_data.index.values,
                                ko_to_pathway_data.index.values))
    ko_abun_data = ko_abun_data.loc[ko]
    ko_to_pathway_data = ko_to_pathway_data.loc[ko]

    ###########################################################################
    # FOR EACH PATHWAY, COUNT THE NUMBER OF NON-ZERO KOS THAT MAP TO IT IN
    # EACH SAMPLE
    ###########################################################################
    ko_binary_data = (ko_abun_data.values > 0).astype(float)
    pathway_counts = np.dot(ko_binary_data.T, ko_to_pathway_data).T

    ###########################################################################
    # CREATE THE KO TO PATHWAY MAPPING
    ###########################################################################
    if args['mapping_method'] != 'naive':
        print("No other method implemented yet...")
        exit()

    ###########################################################################
    # COMPUTE THE PATHWAY ABUNDANCES
    ###########################################################################
    if args['compute_method'] == 'sum':
        pathway_abundance = np.dot(ko_abun_data.T, ko_to_pathway_data).T

    ###########################################################################
    # WRITE OUTPUT
    ###########################################################################

    print("Writing output...")

    if args['transpose_output']:
        path_pd = pd.DataFrame(data=pathway_abundance.T,
                               index=ko_abun_data.columns,
                               columns=ko_to_pathway_data.columns)

    else:
        path_pd = pd.DataFrame(data=pathway_abundance,
                               index=ko_to_pathway_data.columns,
                               columns=ko_abun_data.columns)

    path_pd.index.name = 'Pathway'

    if 'output_file' in args.keys():
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            path_pd.to_csv(args['output_file'], sep='\t')

    elif 'output_pd' in args.keys():
        args['output_pd'] = path_pd

    else:
        sys.exit('Error: No output destination given')

    if 'output_counts_file' in args.keys():
        counts_pd = pd.DataFrame(data=pathway_counts,
                                 index=ko_to_pathway_data.columns,
                                 columns=ko_abun_data.columns)
        counts_pd.index.name = 'Pathway'
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            counts_pd.to_csv(args['output_counts_file'], sep='\t')

    print("Done.")

###############################################################################

if __name__ == "__main__":
    # get options from user
    parser = argparse.ArgumentParser(description='compute the abundance of '
                                                 'pathways from KOs')

    parser.add_argument('-ko', '--ko_abundance', dest='ko_abun_file',
                        help='Input file of ko abundance per sample',
                        default=None)

    parser.add_argument('-ko2path', '--ko_to_pathway',
                        dest='ko_to_pathway_file',
                        help='Input file of mappingfrom ko to pathway',
                        default=None)

    parser.add_argument('-o', '--output', dest='output_file',
                        help='Output file for resulting pathway abundance ('
                             'default: out.tab)', default='out.tab')

    parser.add_argument('-oc', '--output_counts', dest='output_counts_file',
                        help='Output file for number of KOs mapped to each '
                             'pathway (default: counts.tab)',
                        default='counts.tab')

    parser.add_argument('-map', '--mapping_method', dest='mapping_method',
                        help='Method to map KOs to Pathway (default: naive)',
                        default='naive', choices=['naive'])

    parser.add_argument('-compute', '--compute_method', dest='compute_method',
                        help='Method to compute pathway abundance from '
                             'mapped KOs (default: sum)',
                        default='sum', choices=['sum'])

    parser.add_argument('-transpose_ko', '--transpose_ko_abundance',
                        dest='transpose_ko_abundance',
                        help='Transpose the ko abundance matrix given ('
                             'default: False)', action='store_true')

    parser.add_argument('-transpose_output', '--transpose_output',
                        dest='transpose_output',
                        help='Transpose the output pathway abundance matrix '
                             '(default: False)', action='store_true')

    parser.add_argument('-v', '--verbose', dest='verbose',
                        help='Increase verbosity of module (default: false)',
                        action='store_true')

    given_args = parser.parse_args()
    main(vars(given_args))






