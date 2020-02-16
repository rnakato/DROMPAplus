PC ENRICH: Enrichment visualization
-----------------------------------------

For a small genome (e.g., yeast), the sequencing depth is generally enough (> 10 fold).
In such cases, the genome-wide ChIP/Input enrichment distribution is informative because the
technical and biological bias in high throughput sequencing can be minimized.

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

