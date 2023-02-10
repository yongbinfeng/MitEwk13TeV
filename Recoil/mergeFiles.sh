#!/bin/bash
# merge root files based on filenames

INPUTDIR=/eos/uscms/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Zmumu_pT20
OUTPUTDIR=$INPUTDIR/merged

for i in `ls $INPUTDIR | grep root | awk -F'_' '{print $1}' | sort | uniq`
do
    echo $i
    echo "hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*"
    hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*
done
