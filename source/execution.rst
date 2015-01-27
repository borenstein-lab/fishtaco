FiShTaCo API via the command line
===============================
The FiShTaCo module handles all calculations internally.
FiShTaCo offers an interface to the FiShTaCo functionality via the command line and the run_fishtaco.py script.

Usage
------

``run_fishtaco.py input_file [options]``

Required arguments
-------------------

**-ta, --taxa_abundance TAXA_ABUN_FILE**
    Input file of taxa abundance

**-fu, --function_abundance FUNCTION_ABUN_FILE**
    Input file of function abundance

**-c, --class CLASS_FILE**
    Input file of class assignment for the two different
    compared classes


Optional arguments
-------------------

**-h, --help**
    show help message and exit

**-t2f, --taxa_to_function TAXA_TO_FUNCTION_FILE**
    Input file of mapping from taxa to functions

**-op, --output_prefix OUTPUT_PREF**
    Output prefix for result files (default: out)

**-da, --da_results DA_RESULT_FILE**
    Pre-computed DA results from the compute_differential_abundance.py script (default: None)

**-function_da_threshold {Bonf, FDR-0.01, FDR-0.05, FDR-0.1, None}**
    Differential abundance threshold (default: None)

**-max_da MAX_DA_FUNCTIONS_CASES_CONTROLS, --max_da_functions MAX_DA_FUNCTIONS_CASES_CONTROLS**
    Maximum number of differential abundant functions to consider (default: None)

**-assessment {separate_i, permute_all_but_i, permute_only_i, permuted_shapley_orderings}, --taxa_assessment_method {separate_i, permute_all_but_i, permute_only_i, permuted_shapley_orderings}**
    The method used when assessing taxa to compute score (default: permuted_shapley_orderings)

**-score {t_test, mean_diff, median_diff, wilcoxon, log_mean_ratio}, --score_to_compute {t_test, mean_diff, median_diff, wilcoxon, log_mean_ratio}**
    The score to compute for each taxa (default: wilcoxon)

**-max_score MAX_SCORE_CUTOFF, --max_score_cutoff MAX_SCORE_CUTOFF**
    The maximum score cutoff (for example, when dividing by zero) (default: 100)

**-na_rep NA_REP**
    How to represent NAs in the output (default: NA)

**-number_of_permutations NUMBER_OF_PERMUTATIONS**
    number of permutations (default: 100)

**-number_of_shapley_orderings_per_taxa NUMBER_OF_SHAPLEY_ORDERINGS_PER_TAXA**
    number of shapley orderings per taxa (default: 5)

**-use_t2f_as_prior, --use_taxa_to_function_as_prior**
    Learn the taxa copy number of each function, using the given taxa to function file as prior (default: False)

**-residual_mode {as_taxa,remove_residual,as_baseline}**
    How to treat the residual of the functional abundance profile (default: remove_residual)

**-normalization_mode {none,scale_non_permuted,scale_permuted}**
    How to normalize the sample after permuting taxa (default: none)

**-permutation_mode {independent,blocks}**
    How to permute the taxa across samples (default: independent)

**-single_function_filter SINGLE_FUNCTION_FILTER**
    Limit analysis to this single function (default: All)

**-log, --log**
    Write to log file (default: False)

Examples
--------
In the *musicc/examples* directory, the file *simulated_ko_relative_abundance.tab* contains simulated KO abundance measurements of 20 samples described in the
MUSiCC manuscript. Using this file as input for MUSiCC results in the following files:

- simulated_ko_MUSiCC_Normalized.tab (only normalization)
- simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab (normalize and correct using the generic model learned from HMP)
- simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab (normalize and correct learning a new model for each sample)

The commands used were the following (via command line):

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized.tab``

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -c use_generic -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab``

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -c learn_model -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab``
