# Changelog

## 1.7.2 (2020-05-28)
- Modify output filename of PROFILE and MULTICI
- MULTICI: add `--getmaxposi` option to output the max bin value for each site instead of averaged bin value
- Add `--bed12` option to depict BED12 format (mainly for ChromHMM output)

## 1.7.1 (2020-05-26)
- Bug fix in peak-calling (overflow of p_enrich)
- Add script change_chrname_[to|from]_Greek.pl in scripts directory

## 1.7.0 (2020-05-26)
- Adopt CMake to compile

## 1.6.1 (2020-05-25)
- Bug fix: basename(label) of filename in peak-calling
- Switch from <int32_t> to <double> for y-axis bins in pdf

## 1.6.0 (2020-05-21)
- Adopt htslib API for parsing SAM/BAM/CRAM format
- Now available to use on Mac
- Parse2wig+: Add --allchr option that uses all chromosomes to estimate fragment length

## 1.5.0 (2020-05-12)
- Add GENWIG command to generate wig data of ChIP/Input enrichment and p-value distributions
- Add MULTICI command to generate matrix of averaged read density in specified BED sites

## 1.4.1 (2020-04-23)
- Add parsestats4DROMPAplus.pl in scripts directory

## 1.4.0 (2020-03-12)
- Draft the manual page (http://drompaplus.readthedocs.io/)
- Add drompa.heatmap.py for heatmap visualization in otherbins directory
- Add sample scripts and related annotation data in tutorial directory
- Add data/geneannotation directory and add Ensembl gene annotations
- Add parsestats4DROMPAplus.pl and makemappabilitytable.pl in otherbins directory
- Add scripts link in the root directory
- Add symbolic link to cpdf binary file in otherbins/
- Add --shownegative option
- Add --offpdf option
- Set default threshold values for peak calling in each mode
- Add function for highlighting peak regions
- Move the source files from submodule/SSP/common to src
- PROFILE: Bug fix in --stype 2
- Bug fix in ARS visualiztion for yeast
- Bug fix in --showratio 2 (logratio)
- Remove --showars option

## 1.3.2 (2020-02-14)
- Bug fix in --inter for Mango file
- Bug fix in --ideogram

## 1.3.1 (2019-08-22)
- Update submodule SSP

## 1.3.0 (2019-08-18)
- add "data" directory that contains genometable, mptable and ideograms
- add --ideogram option

## 1.2.0 (2019-08-03)
- Parse2wig+: modify --outputzero option from int to bool
- Parse2wig+: add --verbose option. "*.ReadLengthDist.tsv" and "*.ReadCountDist.tsv" are output only when supplying --verbose option.
- DROMPA+: withdraw CI, CG, TR, HEATMAP commands
- DROMPA+: add --width_page and --width_draw options
- DROMPA+: add --chipdrop option to visualize ChIPDrop data
- DROMPA+: allow not to specify --input (-i) option.
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
