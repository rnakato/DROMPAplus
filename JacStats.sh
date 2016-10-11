#!/bin/bash

dir="`pwd`"
/home/rnakato/git/DROMPAplus/parsejaccard.pl
for file in `ls $dir/parse2wigdir+/*.fvp.csv | sed -e 's/.fvp.csv//g'`
do
    echo -en "$file\t"
    /home/rnakato/git/DROMPAplus/parsefvp.pl $file.fvp.csv

#    echo -en "`cut -f2 $file.mpfv.csv | tail -n2 | head -n1`\t"
 #   echo -en "`cut -f3 $file.mpfv.csv | tail -n2 | head -n1`\t"
  #  echo -en "`cut -f4 $file.mpfv.csv | tail -n2 | head -n1`\t"
   # echo -en "`cut -f5 $file.mpfv.csv | tail -n2 | head -n1`\t"
#    echo -en "`cut -f6 $file.mpfv.csv | tail -n2 | head -n1`\t"
#    echo -en "`cut -f7 $file.mpfv.csv | tail -n2 | head -n1`\t"
#    echo -en "`cut -f8 $file.mpfv.csv | tail -n2 | head -n1`\t"
#    echo -en "`cut -f9 $file.mpfv.csv | tail -n2 | head -n1`\t"
#    echo -en "`cut -f10 $file.mpfv.csv | tail -n2 | head -n1`\t"
#    echo -en "`cut -f11 $file.mpfv.csv | tail -n2 | head -n1`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$2} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$3} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$4} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$5} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$6} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$7} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$8} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$9} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$10} END{print m;}'`\t"
#    echo -en "`cat $file.mpfv.csv | awk '{m+=$11} END{print m;}'`\t"
    echo -en "`/home/rnakato/git/DROMPAplus/parsejaccard.pl $file.jaccard.csv`"
#    echo -en "`head -n1 ppout/$file.resultfile | cut -f3,9,10,11`\t"
#    echo -en "`grep "# d = " macs/${file}_peaks.xls | sed -e 's/# d = //g'`\t"                                           
 #   echo -en "`grep "Redundant rate in treatment" macs/${file}_peaks.xls  |  awk -F' ' '{print $6}'`\t"
    #  npeak=`wc -l macs_nomodel/${file}_peaks.xls | awk -F' ' '{print $1}'`
    echo ""
done
