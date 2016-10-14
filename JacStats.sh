#!/bin/bash

dir="`pwd`"
#/home/rnakato/git/DROMPAplus/parsejaccard.pl

echo -en "\tFragment score\tRead score\tlen100\tlen150\tlen500\tlen1000\tlen2000\tlen3000\tlen10000\tlen100000\tlen1000000"
echo -en "\tNSC\tRLSC\tfragment length\tBackground enrichment\tBackground uniformity"
echo -en "\tLen-200\tlen0\tlen100\tlen150\tlen500\tlen1000\tlen1500\tlen5000\tlen100000\tlen200000\tlen300000\tlen400000\tlen500000\tlen600000\tlen700000\tlen800000\tlen900000"
echo -e "\tLen-200\tlen0\tlen100\tlen150\tlen500\tlen1000\tlen1500\tlen5000\tlen100000\tlen200000\tlen300000\tlen400000\tlen500000\tlen600000\tlen700000\tlen800000\tlen900000"

for file in `ls $dir/parse2wigdir+/*.fvp.csv | sed -e 's/.fvp.csv//g'`
do
    echo -en "$file\t"
    parsefvp.pl $file.fvp.csv

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
    echo -en "`parsejaccard.pl -n 3 $file.jaccard.csv`"
    echo -en "`parsejaccard.pl -n 4 $file.jaccard.csv`"
#    echo -en "`head -n1 ppout/$file.resultfile | cut -f3,9,10,11`\t"
#    echo -en "`grep "# d = " macs/${file}_peaks.xls | sed -e 's/# d = //g'`\t"                                           
 #   echo -en "`grep "Redundant rate in treatment" macs/${file}_peaks.xls  |  awk -F' ' '{print $6}'`\t"
    #  npeak=`wc -l macs_nomodel/${file}_peaks.xls | awk -F' ' '{print $1}'`
    echo ""
done
