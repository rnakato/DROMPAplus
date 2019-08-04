# Changelog

## 1.2.0 (2019-08-03)
- Parse2wig+: modify --outputzero option from <int> to <bool>
- Parse2wig+: add --verbose option. "*.ReadLengthDist.tsv" and "*.ReadCountDist.tsv" are output only when supplying --verbose option.
- DROMPA+: withdraw CI, CG, TR, HEATMAP commands
- add `tutorial` directory
- modify help messages

## 1.1.3 (2019-04-15)
- Parse2wig+: Accept CRAM format for the input file
- src/SSP/script/makegenometable.pl: accept scaffold in Ensembl genome

## 1.1.2 (2018-08-29)
- bug fix for pdf page count

## 1.1.1 (2018-08-24)
- modify PROFILE command for bedsites

## 1.1.0 (2018-07-22)
- remove ZINB estimation for peak-calling
- change bindist.csv to readcountdist.tsv
- change to use the longest chromosome only for estimating fragment length in default

## 1.0.0
- First commit
