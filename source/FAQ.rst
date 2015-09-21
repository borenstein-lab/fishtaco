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

