#!/bin/bash

dir=parse2wigdir+
/home/rnakato/git/DROMPAplus/parsejaccard.pl 
for file in `ls $dir/*.jaccard.csv | sed -e 's/'$dir'\///g' -e 's/.jaccard.csv//g'`
do
    echo -en "$file\t"
    echo -en "`grep "Fragment score" $dir/$file.fvp.csv |cut -f2`\t"
    echo -en "`grep "Read score" $dir/$file.fvp.csv |cut -f2`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$2} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$3} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$4} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$5} END{print m;}'`\t"
    echo -en "`cat $dir/$file.mpfv.csv | awk '{m+=$6} END{print m;}'`\t"
    /home/rnakato/git/DROMPAplus/parsejaccard.pl $dir/$file.jaccard.csv
done
