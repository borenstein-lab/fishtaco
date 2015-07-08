FishTaco Visualization
======================

.. index:: Visualization

The functional shift decomposition figures shown in the FishTaco paper are produced in R (http://www.r-project.org/) by using the FishTacoPlot package.

FishTacoPlot package
--------------------

Plotting the results obtained from applying FishTaco to you data is done via the FishTacoPlot package.

We used R for the plotting in order to take full advantage of the capabilities of the `ggplot2 package <http://ggplot2.org/>`_.

We recommend downloading and installing `R Studio <http://www.rstudio.com/>`_ before plotting, but it is not required.

The following packages need to be pre-installed in your R environment before using the FishTacoPlot package:

* package-scales {scales}
* ggplot2 {ggplot2}

In order to use the FishTacoPlot package:

1. Download the FishTacoPlot package from `GitHub <https://github.com/omanor/fishtaco-plot/archive/1.0.0.tar.gz>`_.

2. Install the package in your R terminal with the command: ``install.packages(<path_to_package>, repos = NULL, type="source")``

3. Load the packages to your workspace with the command: ``require(FishTacoPlot); require(ggplot2); require(scales)``

4. Plot your results by using the function *MultiFunctionTaxaContributionPlots*.

Examples
--------

In the *<path_to_package>/examples* directory, you can find results obtained by applying FishTaco to a healthy cohort of different body sites (HMP),
and a type 2 diabetes cohort (T2D).

Metagenomic Comparative Analysis of Different Body Sites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example of a FishTaco plot showing the species-level contributions to the functional enrichment of 3 pathways in tongue samples compared to
cheek samples:

.. figure:: FishTaco_HMP.png
    :width: 750px
    :align: center
    :height: 500px
    :alt: alternate text
    :figclass: align-center


Generate this example by using the following command:

.. code:: python

    MultiFunctionTaxaContributionPlots(input_dir="<path_to_package>/examples", input_prefix="HMP_fishtaco",
    input_taxa_taxonomy="<path_to_package>/examples/HMP_TAXONOMY.tab", sort_by="list", plot_type="bars",
    input_function_filter_list=c("ko00020", "ko00540","ko02040"), add_predicted_da_markers=TRUE, add_original_da_markers=TRUE)


Metagenomic Comparative Analysis of a Disease Cohort
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example of a FishTaco plot showing the genus-level contributions to the functional enrichment of 3 modules in type 2 diabetes patients compared to
healthy controls:

.. figure:: FishTaco_T2D.png
    :width: 750px
    :align: center
    :height: 500px
    :alt: alternate text
    :figclass: align-center


Generate this example by using the following command:

.. code:: python

    MultiFunctionTaxaContributionPlots(input_dir="<path_to_package>/examples", input_prefix="T2D_fishtaco",
    input_taxa_taxonomy="<path_to_package>/examples/T2D_TAXONOMY.tab", sort_by="list", plot_type="bars",
    input_function_filter_list=c("ko00020", "ko00540","ko02040"), add_predicted_da_markers=TRUE)
