Running FishTaco
================

.. index:: Execution

The FishTaco python module handles all calculations internally.
FishTaco offers an interface to the FishTaco functionality via the command line and the run_fishtaco.py script.

Usage
-----

``run_fishtaco.py -ta TAXA_ABUN_FILE -fu FUNCTION_ABUN_FILE -l LABELS_FILE [options]``

Required arguments
------------------

``-ta, --taxa_abundance TAXA_ABUN_FILE``
    Input file of taxonomic abundance profiles (`format <fishtaco_file_formats.html#taxa-abundance-file>`_)

``-fu, --function_abundance FUNCTION_ABUN_FILE``
    Input file of function abundance (`format <fishtaco_file_formats.html#function-abundance-file>`_)

``-l, --labels LABELS_FILE``
    Input file of label assignment for the two sample sets being compared (`format <fishtaco_file_formats.html#sample-sets-labels-file>`_)


Recommended arguments
---------------------

``-gc, --genomic_content_of_taxa GENOMIC_CONTENT_FILE``
    Input file of genomic content of each taxa (`format <fishtaco_file_formats.html#genomic-content-file>`_)

``-control_label LABEL``
    Define control label (default: 0)

``-case_label LABEL``
    Define case label (default: 1)

``-op, --output_prefix OUTPUT_PREF``
    Output prefix for result files (default: fishtaco_out)

``-map_function_level [pathway,module,none,custom]``
    Map KOs to pathways, modules, none, or custom (default: pathway)


Advanced usage arguments
------------------------

``-map_function_file FUNC_LEVEL_MAP_FILE``
    Mapping file from KOs to pathways, modules, or custom (default: use KEGG database downloaded 07/15/2013)

``-mult, --multiple_hypothesis_correction [Bonf,FDR-0.01,FDR-0.05,FDR-0.1,none]``
    Multiple hypothesis correction for functional enrichment (default: FDR-0.05)

``-max_func, --maximum_functions_to_analyze MAX_FUNCTIONS``
    Maximum number of enriched functions to consider (default: All)

``-assessment, --taxa_assessment_method [permute_all_but_i,permuted_shapley_orderings]``
    The method used when assessing taxa to compute score (default: permuted_shapley_orderings)

``-score, --score_to_compute [t_test,mean_diff,median_diff,wilcoxon,log_mean_ratio]``
    The enrichment score to compute for each function (default: wilcoxon)

``-max_score, --max_score_cutoff MAX_SCORE_CUTOFF``
    The maximum score cutoff (for example, when dividing by zero) (default: 100)

``-na_rep NA_REP``
    How to represent NAs in the output (default: NA)

``-number_of_permutations NUMBER_OF_PERMUTATIONS``
    number of permutations (default: 100)

``-number_of_shapley_orderings_per_taxa NUMBER_OF_SHAPLEY_ORDERINGS_PER_TAXA``
    number of shapley orderings per taxa (default: 5)

``-use_gc_as_prior, --use_genomic_content_of_taxa_as_prior``
    Learn the taxa copy number of each function, using the given genomic content data as prior (default: False)

``-residual_mode [as_taxa,remove_residual,as_baseline]``
    How to treat the residual of the functional abundance profile (default: remove_residual)

``-normalization_mode [none,scale_non_permuted,scale_permuted]``
    How to normalize the sample after permuting taxa (default: scale_permuted)

``-permutation_mode [independent,blocks]``
    How to permute the taxa across samples (default: blocks)

``-en, --enrichment_results DA_RESULT_FILE``
    Pre-computed functional enrichment results from the compute_differential_abundance.py script (default: None)

``-single_function_filter SINGLE_FUNCTION_FILTER``
    Limit analysis only to this single function (default: None)

``-h, --help``
    show help message and exit

``-functional_profile_already_corrected_with_musicc``
    Indicates that the functional profile has been already corrected with MUSiCC prior to running FishTaco (default: False)

``-log, --log``
    Write to log file (default: False)


FishTaco Output Files
---------------------

Main output files
^^^^^^^^^^^^^^^^^
``fishtaco_out_main_output_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the taxon-level decomposition of shift scores for the differentially abundant functions. (`format <fishtaco_file_formats.html#fishtaco-output-file>`_)

Supporting stats output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``fishtaco_out_STAT_taxa_contributions_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the final taxon-level contribution score for every differentially abundant(shifted) function in the input data, as calculated by FishTaco

``fishtaco_out_STAT_DA_function_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains statistics regarding the differential abundance for each function in the input file

``fishtaco_out_STAT_DA_taxa_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains statistics regarding the differential abundance for each taxa in the input file

``fishtaco_out_STAT_mean_stat_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the mean taxon-level contribution score for every differentially abundant(shifted) function in the input data (in default settings, this is equal to the final score)

``fishtaco_out_STAT_median_stat_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the median taxon-level contribution score for every differentially abundant(shifted) function in the input data

``fishtaco_out_STAT_std_stat_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the standard deviation of taxon-level contribution score for every differentially abundant(shifted) function in the input data

``fishtaco_out_STAT_original_value_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the metagenome-based shift statistics value for each function in the input file

``fishtaco_out_STAT_predicted_DA_value_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the taxa-based shift statistics value for each function in the input file

``fishtaco_out_STAT_predicted_function_abundance_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the taxa-based abundance profile for each function in each sample

``fishtaco_out_STAT_predicted_function_agreement_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains various statistics regarding the agreement between the metagenome- and taxa-based abundance profiles for each function

``fishtaco_out_STAT_residual_function_abundance_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the residual between the metagenome- and taxa-based abundance profiles for each function (in "remove-residual" mode the residual is equal to zero)

``fishtaco_out_STAT_shapley_orderings_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the random Shapley orderings used in the run (for "permuted_shapley_orderings" mode)

``fishtaco_out_STAT_taxa_learned_copy_num_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the inferred copy numbers of each function in each taxon (for FishTaco with prior-based or *de novo* inference)

``fishtaco_out_STAT_taxa_learning_rsqr_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains various statistics regarding the agreement between the metagenome- and taxa-based abundance profiles for each function (on test data)

``fishtaco_out_STAT_run_log_SCORE_wilcoxon_ASSESSMENT_permuted_shapley_orderings.tab``
    contains the running log of FishTaco

Examples
--------
The *fishtaco/examples* directory contains the following files:

- the file *METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab* contains scaled abundance measurements of 10 species in 213 samples from the HMP dataset
- the file *WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab* contains MUSiCC-corrected abundance values for the K00001 orthology group in the same samples
- the file *METAPHLAN_taxa_vs_KO_only_K00001.tab* contains the copy numbers of the K00001 orthology group in the 10 species as above
- the file *SAMPLE_vs_CLASS.tab* contains class labels from the same samples (control vs. case)

Using these files as input for FishTaco results in the following output files (found in the *fishtaco/examples/output* directory):

FishTaco with no inference
^^^^^^^^^^^^^^^^^^^^^^^^^^

Running FishTaco with no inference generates the output files found in output/fishtaco_out_no_inf_STAT_*

.. code:: python

    run_fishtaco.py -ta fishtaco/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab
    -fu fishtaco/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab
    -t2f fishtaco/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab
    -c fishtaco/examples/SAMPLE_vs_CLASS.tab -op fishtaco_out_no_inf -max_da 1
    -assessment permuted_shapley_orderings -score wilcoxon -na_rep 0
    -number_of_shapley_orderings_per_taxa 3 -residual_mode remove_residual -log
    -normalization_mode scale_permuted -permutation_mode blocks -number_of_permutations 5

FishTaco with prior-based inference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running FishTaco with prior-based inference generates the output files found in output/fishtaco_out_prior_based_inf_STAT_*


.. code:: python

    run_fishtaco.py -op fishtaco_out_no_inf -max_da 1
    -ta fishtaco/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab
    -fu fishtaco/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab
    -c fishtaco/examples/SAMPLE_vs_CLASS.tab
    -t2f fishtaco/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab
    -assessment permuted_shapley_orderings -score wilcoxon
    -na_rep 0 -number_of_shapley_orderings_per_taxa 3 -residual_mode remove_residual
    -normalization_mode scale_permuted -permutation_mode blocks -number_of_permutations 5
    -use_t2f_as_prior -log

FishTaco with de novo inference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running FishTaco with *de novo* inference generates the output files found in output/fishtaco_out_de_novo_inf_STAT_*


.. code:: python

    run_fishtaco.py -op fishtaco_out_no_inf -max_da 1
    -ta fishtaco/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab
    -fu fishtaco/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab
    -c fishtaco/examples/SAMPLE_vs_CLASS.tab -assessment permuted_shapley_orderings
    -score wilcoxon -na_rep 0 -number_of_shapley_orderings_per_taxa 3
    -residual_mode remove_residual -normalization_mode scale_permuted
    -permutation_mode blocks -number_of_permutations 5 -log































