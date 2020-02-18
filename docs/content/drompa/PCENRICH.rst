PC ENRICH: Enrichment visualization
-----------------------------------------

For a small genome (e.g., yeast), the sequencing depth is generally enough (> 10 fold).
In such cases, the genome-wide ChIP/Input enrichment distribution is informative because the
technical and biological bias in high throughput sequencing can be minimized.


Download the data
+++++++++++++++++++++++++++++++

Download the sample data (CRAM-format map files) from GoogleDrive::

    https://drive.google.com/open?id=1f5H-umPgjzmDYLVHBdlIMWQXtt91S6Fc
    https://drive.google.com/open?id=1f991hi-V9ITAWeUzr3AkqyRY8CbFTBXO
    https://drive.google.com/open?id=1fDUgeo8IvhI-fikQv_PDCiDmNc9cBfWT
    https://drive.google.com/open?id=1fGjY5nlhjePk_0TXKC5bLOzmdrer7IHh
    https://drive.google.com/open?id=1fJh-f3CDNc7wAkuN79rndjptVq_qYU6k
    https://drive.google.com/open?id=1fNRt1uvA1CQrIfb9NSJq3hYL9-GGJnXv


Make enrichment distribution
++++++++++++++++++++++++++++++++++++++++++


To make a PDF file of the enrichment distribution for S. cerevisiae, type::

  $ drompa+ PC_ENRICH -p drompa4 $s1 $s2 $s3 $s4 $s5 $s6 \
  $ --gt genometable.txt -g SGD_features.tab --gftype 3 --ars ARS-oriDB.txt \
  $ --scale_ratio 4 --ls 100

``--ars`` option is specificied to visualize DNA replication origin (ARS) available for *S. cerevisiae* and *S. pombe*. The annotation data can be obtained from `OriDB <http://cerevisiae.oridb.org/>`_.

Supply ``--showratio 2`` to show logratio::

  $ drompa+ PC_ENRICH -p drompa4 $s1 $s2 $s3 $s4 $s5 $s6 \
  $ --gt genometable.txt -g SGD_features.tab --gftype 3 --ars ARS-oriDB.txt \
  $ --scale_ratio 4 --ls 100 --showratio 2



If you want to check the enrichment precisely, it is good to adjust the y-axis as follows::

  $ drompa+ PC_ENRICH -p drompa5 $s1 $s2 $s3 $s4 $s5 $s6 \
  $ --gt genometable.txt -g SGD_features.tab --gftype 3 -ars ARS-oriDB.txt \
  $ --scale_ratio 4 --ls 100 -bn 5 -ystep 10

