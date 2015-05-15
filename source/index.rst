.. FishTaco documentation master file, created by
   sphinx-quickstart on Mon Jan 26 15:16:56 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FishTaco: Functional Shift Taxonomic Contributors
=================================================
.. index:: Home

FishTaco is a metagenomic computational framework, aiming to identify the taxa that are driving the functional shifts
we observe in microbiomes of different individuals or disease states.

.. figure:: FishTaco_pipeline.png
    :width: 750px
    :align: center
    :height: 500px
    :alt: alternate text
    :figclass: align-center

    The FishTaco framework. (A) Deconvolving the metagenome-based functional abundance profile into a taxa-based
    functional profile, based on the taxonomic profile and the genomic content of taxa (left).
    The genomic content is an optional input as FishTaco can infer the genomic content from the taxonomic and
    functional profiles. The metagenome- and taxa-based functional shifts are calculated by comparative functional
    analysis (right). (B) Decomposing the functional shifts into taxon-level contributions. By using a permutation-based
    approach, the context-dependent contribution of a single taxon to the functional shift is quantified.
    To account for statistical interactions between taxa and to linearize the contribution profile, a Shapley value analysis is used.
    (C) The FishTaco taxon-level contribution profile. Each function’s functional shift score is decomposed into taxon-level
    contributions (left), distinguishing between 4 modes of contribution (right).


.. toctree::
   :maxdepth: 2

   installation
   execution
   visualization


==================
Citing Information
==================

If you use the FishTaco software, please cite the following paper:

**Linking taxonomic and functional metagenomic profiles to identify species driving functional shifts in the human microbiome.**
Ohad Manor and Elhanan Borenstein. *In preparation*

==================
Question forum
==================
For FishTaco announcements and questions, including notification of new releases, you can visit the `FishTaco users forum <https://groups.google.com/forum/#!forum/fishtaco-users>`_.

