Frequently asked questions
==========================
.. index:: FAQ





Can I run FishTaco with PICRUSt-derived metagenomic functional profile?
-----------------------------------------------------------------------
Yes, FishTaco can be also used for 16S and PICRUSt-predicted metagenomes. The user needs to supply the taxonomic abundance data
(output of 16S profiling), combined with the functional profile (PICRUSt-predicted metagenome), and optionally,
the PICRUSt pre-calculated files for genomic content
(available in the (`PICRUSt website <http://picrust.github.io/picrust/picrust_precalculated_files.html#id1>`_)


What is the running time for FishTaco?
--------------------------------------
The running time of FishTaco varies greatly depending on two factors: the number of functions examined, and the number of taxa in each sample.
Since FishTaco uses a permutation-based approach to quantify the contribution of each taxon to every functional shift, the running time is linearly
dependent on the number of taxa. When running on a single CPU process, decomposing the functional shift associated with a
single function in a dataset with ~100 taxa took 2-3 hours on average, while running on a dataset with ~50 taxa took 1-2 hours on average. Therefore,
if you are running FishTaco on a dataset with a large number of taxa and functions, please read the answer to the question about making FishTaco
run faster.


How can I make FishTaco run faster?
-----------------------------------
In order to get a faster decomposition of the functional shift to the taxon-level contributions, you can try the following:

Reducing the number of taxa
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Since FishTaco running time depends on the number of taxa, you can remove taxa that are rare in terms if their average abundance across all
samples (as they are less likely to drive a significant functional change) or in terms of the number of samples they are present in
(for the same reason).

Reducing the number of functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Since FishTaco running time depends on the number of functions, you can run FishTaco on a reduced number of functions, by either defining a more
strict cutoff for significance of functional shift (since FishTaco only decomposes significantly shifted functions, use the flag *-mult_hyp*),
or running in the default pathway-level mode as opposed to KO or Module modes. In addition, since FishTaco can be ran on a single function at a time
(using the flag *-single_function_filter*), and the functions are independent, advanced users can run FishTaco in parallel for all their functions of
interest and merge the results.

Running in single-taxa mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^
By default, FishTaco runs in a multi-taxa mode, where the contributions of all taxa are quantified simultaneously. However, using the flag
*-assessment single_taxa* will result in much faster (x100) but a less-accurate decomposition of the functional shifts
(this for example can be used to select interesting functions for multi-taxa analysis).

