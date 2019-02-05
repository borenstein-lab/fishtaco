Running FishTaco
================

.. index:: Execution

The FishTaco python module handles all calculations internally.
FishTaco offers an interface to the FishTaco functionality via the command line and the run_fishtaco.py script.

Usage
-----

FishTaco can be used in two alternative modes, depending on the availability of **genomic information** for each taxon. Specifically,
if such data is available (e.g., through reference genomes), FishTaco can be used with the *-gc* flag. However, FishTaco can also infer this
data by using the *-inf* flag. If you are using 16S data coupled with PICRUSt,
please read :ref:`picrust-info`.

**Running FishTaco with genomic content data:**

``run_fishtaco.py -ta TAXA_ABUN_FILE -fu FUNCTION_ABUN_FILE -l LABELS_FILE -gc GENOMIC_CONTENT_FILE [options]``

**Running FishTaco with genomic content inference:**

``run_fishtaco.py -ta TAXA_ABUN_FILE -fu FUNCTION_ABUN_FILE -l LABELS_FILE -inf [options]``


Required arguments
------------------

``-ta, --taxa_abundance TAXA_ABUN_FILE``
    Input file of taxonomic abundance profiles (`format <fishtaco_file_formats.html#taxa-abundance-file>`_)

``-l, --labels LABELS_FILE``
    Input file of label assignment for the two sample sets being compared (`format <fishtaco_file_formats.html#sample-sets-labels-file>`_)


Optional arguments
------------------

``-fu, --function_abundance FUNCTION_ABUN_FILE``
    Input file of function abundance (`format <fishtaco_file_formats.html#function-abundance-file>`_)

``-gc, --genomic_content_of_taxa GENOMIC_CONTENT_FILE``
    Input file of genomic content of each taxa (`format <fishtaco_file_formats.html#genomic-content-file>`_)

``-inf, --perform_inference_of_genomic_content``
    Defines if genome content is inferred (either de-novo or prior-based if genomic content is also given, default: FALSE)

``-label_to_find_enrichment_in``
    Define sample set label to find enrichment in (default: 1)

``-label_to_find_enrichment_against``
    Define sample set label to find enrichment against (default: 0)

``-op, --output_prefix OUTPUT_PREF``
    Output prefix for result files (default: fishtaco_out)

``-map_function_level {pathway, module, none, custom}``
    Map KOs to pathways, modules, none, or custom (default: pathway)

``-assessment, --taxa_assessment_method {single_taxa, multi_taxa}``
    The method used when assessing taxa to compute individual contributions. The running time of *single_taxa* will
    be significantly lower than *multi_taxa*, but less accurate (see manuscript for details) (default: multi_taxa)

Advanced usage arguments
------------------------

``-map_function_file FUNC_LEVEL_MAP_FILE``
    Mapping file from KOs to pathways, modules, or custom (default: use internal KEGG database downloaded 07/15/2013)

``-perform_inference_on_ko_level``
    Indicates to perform the inference on the KO level (default: use the mapped functional level, e.g., pathway)

``-mult_hyp, --multiple_hypothesis_correction {Bonf, FDR-0.01, FDR-0.05, FDR-0.1, none}``
    Multiple hypothesis correction for functional enrichment (default: FDR-0.05)

``-max_func, --maximum_functions_to_analyze MAX_FUNCTIONS``
    Maximum number of enriched functions to consider (default: All)

``-score, --score_to_compute {t_test, mean_diff, median_diff, wilcoxon, log_mean_ratio}``
    The enrichment score to compute for each function (default: wilcoxon)

``-max_score, --max_score_cutoff MAX_SCORE_CUTOFF``
    The maximum score cutoff (for example, when dividing by zero) (default: 100)

``-na_rep NA_REP``
    How to represent NAs in the output (default: NA)

``-number_of_permutations NUMBER_OF_PERMUTATIONS``
    number of permutations (default: 100)

``-number_of_shapley_orderings_per_taxa NUMBER_OF_SHAPLEY_ORDERINGS_PER_TAXA``
    number of shapley orderings per taxa (default: 5)

``-en, --enrichment_results DA_RESULT_FILE``
    Pre-computed functional enrichment results from the compute_differential_abundance.py script (default: None)

``-single_function_filter SINGLE_FUNCTION_FILTER``
    Limit analysis only to this single function (default: None)

``-multi_function_filter_list MULTI_FUNCTION_FILTER``
    Limit analysis only to these comma-separated functions (default: None)

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
    contains the taxon-level decomposition of shift scores for the differentially abundant functions. (`format <fishtaco_file_formats.html#fishtaco-output-file-main-output>`_)

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
The *fishtaco/examples* directory contains the following:

- the file *METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab* contains scaled abundance measurements of 10 species in 213 samples from the HMP dataset
- the file *WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab* contains MUSiCC-corrected abundance values for the K00001 orthology group in the same samples
- the file *METAPHLAN_taxa_vs_KO_only_K00001.tab* contains the copy numbers of the K00001 orthology group in the 10 species as above
- the file *SAMPLE_vs_CLASS.tab* contains class labels from the same samples (control vs. case)

Using these files as input for FishTaco results in the following output files (found in the *fishtaco/examples/output* directory):

Note: If you installed the FishTaco package using *pip*, the *examples* directory is located in your python packages directory, e.g.,
*lib/python3.3/site-packages*

FishTaco with no inference
^^^^^^^^^^^^^^^^^^^^^^^^^^

Running FishTaco with no inference generates the output files found in fishtaco/examples/output/fishtaco_out_no_inf_STAT_*

.. code:: python

    run_fishtaco.py -ta fishtaco/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab -fu fishtaco/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab
    -l fishtaco/examples/SAMPLE_vs_CLASS.tab -gc fishtaco/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab -op fishtaco_out_no_inf
    -map_function_level none -functional_profile_already_corrected_with_musicc -assessment single_taxa -log

FishTaco with prior-based inference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running FishTaco with prior-based inference generates the output files found in fishtaco/examples/output/fishtaco_out_prior_based_inf_STAT_*


.. code:: python

    run_fishtaco.py -ta fishtaco/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab -fu fishtaco/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab
    -l fishtaco/examples/SAMPLE_vs_CLASS.tab -gc fishtaco/examples/METAPHLAN_taxa_vs_KO_only_K00001.tab -op fishtaco_out_prior_based_inf
    -map_function_level none -functional_profile_already_corrected_with_musicc -inf -assessment single_taxa -log

FishTaco with de novo inference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running FishTaco with *de novo* inference generates the output files found in fishtaco/examples/output/fishtaco_out_de_novo_inf_STAT_*


.. code:: python

    run_fishtaco.py -ta fishtaco/examples/METAPHLAN_taxa_vs_SAMPLE_for_K00001.tab -fu fishtaco/examples/WGS_KO_vs_SAMPLE_MUSiCC_only_K00001.tab
    -l fishtaco/examples/SAMPLE_vs_CLASS.tab -op fishtaco_out_de_novo_inf -map_function_level none -functional_profile_already_corrected_with_musicc
    -inf -assessment single_taxa -log


.. raw:: html

    <br>

.. raw:: html

    <br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>



























