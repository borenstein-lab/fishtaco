=======
HISTORY
=======

====================
1.1.5 (24 May, 2021)
====================
* fixed type in read_csv file separator

========================
1.1.4 (5 November, 2020)
========================
* fixed bug that could cause OTU IDs to be read in as integers (caused problems when setting the dataframe index). Due to an underlying pandas bug, could be reverted if pandas is updated

=======================
1.1.3 (23 August, 2019)
=======================
* fixed bug in specifying FDR correction
* changed DA function statistics output file to include Bonferroni- and FDR-corrected significance values, rather than  just whether they passed the supplied alpha threshold
* added tests for FDR correction filtering and changing the alpha cutoff for differentially abundant functions

=======================
1.1.2 (22 August, 2019)
=======================
* fixed deprecated scikit-learn imports
* fixed swapped control and case sample numbers reported while running (this bug was only a display error, it had no effect on the results)
* added option to specify the number of folds for cross validation during genomic content inference
* added informative error message for when there are too few control and/or case samples for genomic content inference (given the number of folds)
* added option to specify the significance cutoff for multiple hypothesis corrections

========================
1.1.1 (20 October, 2016)
========================
* improved unit testing to include each individual function, as well as comparing the calculated FishTaco output to a pre-computed example output
* removed debugging statements from code
* added documentation to all functions
* resolved several runtime warnings for edge cases

=========================
1.1.0 (18 February, 2016)
=========================
* added the option to run FishTaco with no functional profile input, predicting it from the taxonomic profiles and genomic content inputs
* added a unit test for the option to run with no functional profiles input
* added the option to run FishTaco on only a subset of functions by supplying a list file as input

=====================
1.0.5 (15 July, 2015)
=====================
* Initial release

