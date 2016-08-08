#!/bin/bash

dir=parse2wigdir+
/home/rnakato/git/DROMPAplus/parsejaccard.pl
for file in `ls $dir/*.mpfv.csv | sed -e 's/'$dir'\///g' -e 's/.mpfv.csv//g'`
do
    echo -en "$file\t"
    /home/rnakato/git/DROMPAplus/parsefvp.pl $dir/$file.fvp.csv

    echo -en "`cut -f2 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f3 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f4 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f5 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f6 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f7 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f8 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f9 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f10 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cut -f11 $dir/$file.mpfv.csv | tail -n2 | head -n1`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$2} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$3} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$4} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$5} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$6} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$7} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$8} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$9} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$10} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$11} END{print m;}'`\t"
    /home/rnakato/git/DROMPAplus/parsejaccard.pl $dir/$file.jaccard.csv
done
