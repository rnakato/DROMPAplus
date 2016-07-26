#!/bin/bash

/home/rnakato/git/DROMPAplus/parsejaccard.pl 
for file in `ls parse2wigdir/*.jaccard.csv | sed -e 's/parse2wigdir\///g' -e 's/.jaccard.csv//g'`
do
    echo -en "$file\t"
    /home/rnakato/git/DROMPAplus/parsejaccard.pl parse2wigdir/$file.jaccard.csv
done
