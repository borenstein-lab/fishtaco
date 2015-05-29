FishTaco file formats
=====================


Format of taxa abundance file
-----------------------------

A tab-delimited file (with headers), where each row is a taxon and each column is a sample,
and each cell represents the relative abundance of the taxon in the sample:

========  ======== ======== ========
  Taxa    Sample1  Sample2  Sample3
========  ======== ======== ========
*Taxon1*    0.2    0.5      0.1
*Taxon2*    0.3    0.1      0.7
*Taxon3*    0.3    0.2      0.1
*Taxon4*    0.2    0.2      0.1
========  ======== ======== ========

Format of function abundance file
---------------------------------

A tab-delimited file (with headers), where each row is a function (e.g., KO) and each column is a sample,
and each cell represents the relative abundance of the function in the sample:

========  ======== ======== ========
Function  Sample1  Sample2  Sample3
========  ======== ======== ========
*K00001*    0.1    0.1      0.1
*K00002*    0.3    0.1      0.1
*K00003*    0.1    0.1      0.2
*K00004*    0.2    0.1      0.2
*K00005*    0.2    0.1      0.2
*K00006*    0.1    0.5      0.2
========  ======== ======== ========

Format of genome content file
-----------------------------

A tab-delimited file (with headers), where each row is a taxon and each column is a function,
and each cell represents the copy number of the function in the genome of the taxon:

========  ======== ======== ======== ======== ======== ========
Taxa       K00001   K00002  K00003    K00004    K00005  K00006
========  ======== ======== ======== ======== ======== ========
*Taxon1*    1       2           0        1       2          0
*Taxon1*    0       3           0        1       0          0
*Taxon1*    2       1           2        2       3          1
*Taxon1*    0       0           2        0       3          2
========  ======== ======== ======== ======== ======== ========

