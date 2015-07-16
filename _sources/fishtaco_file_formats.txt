FishTaco file formats
=====================

Taxa abundance file
-------------------

A tab-delimited text file (with headers), where each row is a taxon and each column is a sample,
and each cell represents the relative abundance of the taxon in the sample:

.. raw:: html

    <br>

.. rst-class:: center-align

==========  ======== ======== ========
  Taxa      Sample1  Sample2  Sample3
==========  ======== ======== ========
**Taxon1**    0.2    0.5      0.1
**Taxon2**    0.3    0.1      0.7
**Taxon3**    0.3    0.2      0.1
**Taxon4**    0.2    0.2      0.1
==========  ======== ======== ========

.. raw:: html

    <br><br><br>


Function abundance file
-----------------------

A tab-delimited text file (with headers), where each row is a function (e.g., KO) and each column is a sample,
and each cell represents the relative abundance of the function in the sample:

.. raw:: html

    <br>

.. rst-class:: center-align

==========  ======== ======== ========
Function     Sample1  Sample2  Sample3
==========  ======== ======== ========
**K00001**    0.1    0.1      0.1
**K00002**    0.3    0.1      0.1
**K00003**    0.1    0.1      0.2
**K00004**    0.2    0.1      0.2
**K00005**    0.2    0.1      0.2
**K00006**    0.1    0.5      0.2
==========  ======== ======== ========

.. raw:: html

    <br><br><br>

Sample sets labels file
-----------------------

A tab-delimited text file (with headers), with two columns, where each row is a sample the value is the sample label:

.. raw:: html

    <br>

.. rst-class:: center-align

===========  ========
Sample        Label
===========  ========
**Sample1**     1
**Sample2**     0
**Sample3**     1
===========  ========

.. raw:: html

    <br><br><br>

Genomic content file
--------------------

A tab-delimited text file (with headers), where each row is a taxon and each column is a function,
and each cell represents the copy number of the function in the genome of the taxon:

.. raw:: html

    <br>

.. rst-class:: center-align

==========  ======== ======== ======== ======== ======== ========
Taxa         K00001   K00002  K00003    K00004    K00005  K00006
==========  ======== ======== ======== ======== ======== ========
**Taxon1**    1       2           0        1       2          0
**Taxon2**    0       3           0        1       0          0
**Taxon3**    2       1           2        2       3          1
**Taxon4**    0       0           2        0       3          2
==========  ======== ======== ======== ======== ======== ========

.. raw:: html

    <br><br><br>

FishTaco output file (*_main_output_*)
--------------------------------------

A tab-delimited text file (with headers), where each row is a taxon and each column is a function,
and each cell represents the mode and value of contribution to the observed functional shift:

.. raw:: html

    <br>

.. rst-class:: center-align

==========  ======== ======== ======== ======== ======== ========
Taxa         K00001   K00002  K00003    K00004   K00005   K00006
==========  ======== ======== ======== ======== ======== ========
**Taxon1**  a:0.038  a:0.4    c:0.01   c:0.34   b:-0.004 d:-0.01
**Taxon2**  c:0.293  b:-0.03  d:-0.002 c:0.003  a:0.0631 a:0.2139
**Taxon3**  b:-0.008 b:-0.015 a:0.03   c:0.0008 d:-0.009 d:-0.011
**Taxon4**  a:0.048  d:-0.061 b:-0.004 a:0.019  c:4.996  c:0.045
==========  ======== ======== ======== ======== ======== ========

.. raw:: html

    <br>

**a**: case-associated taxon driving enrichment of function in cases;
**b**: case-associated taxon attenuating enrichment of function in cases;
**c**: control-associated taxon driving enrichment of function in cases;
**d**: control-associated taxon attenuating enrichment of function in cases;

.. raw:: html

    <br><br><br>

FishTacoPlot taxonomy file (for visualization)
----------------------------------------------

In order to plot FishTaco's decomposition of functional shifts into taxon-level contributions, a taxonomy file
linking the various taxa to their phylogenetic assignment is required. The taxonomy file is a tab-delimited
text file with **no** headers, where each row represents one taxon, followed by its assignment to the following
phylogenetic levels: 'kingdom','phylum','class','order','family','genus','species', all prefixed by the first letter
of the level followed by 2 underscores.

.. raw:: html

    <br>

.. rst-class:: center-align

==========================  =========== ================= ======================= =================== ===================== ================ =================
Abiotrophia defectiva	    k__Bacteria	p__Firmicutes	  c__Bacilli	          o__Lactobacillales  f__Aerococcaceae	    g__Abiotrophia	 s__defectiva
Acaryochloris marina	    k__Bacteria	p__Cyanobacteria  c__Cyanophyceae	      o__Synechococcales  f__Acaryochloridaceae	g__Acaryochloris s__marina
Acetivibrio cellulolyticus	k__Bacteria	p__Firmicutes	  c__Clostridia	          o__Clostridiales	  f__Clostridiaceae	    g__Acetivibrio	 s__cellulolyticus
Acetobacter pasteurianus	k__Bacteria	p__Proteobacteria c__Alphaproteobacteria  o__Rhodospirillales f__Acetobacteraceae	g__Acetobacter	 s__pasteurianus
Acetohalobium arabaticum	k__Bacteria	p__Firmicutes	  c__Clostridia	          o__Halanaerobiales  f__Halobacteroidaceae	g__Acetohalobium s__arabaticum
==========================  =========== ================= ======================= =================== ===================== ================ =================

+------------------------+------------+----------+----------+
| Header row, column 1   | Header 2   | Header 3 | Header 4 |
| (header rows optional) |            |          |          |
+------------------------+------------+----------+----------+
| body row 1, column 1   | column 2   | column 3 | column 4 |
+------------------------+------------+----------+----------+
| body row 2             | ...        | ...      |          |
+------------------------+------------+----------+----------+

.. raw:: html

    <br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>




































