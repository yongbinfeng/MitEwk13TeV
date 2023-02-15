#!/bin/bash
# merge root files based on filenames

doZmm=True
doZee=True
doWm=True
doWe=True

if [ $doZmm = True ]; then
    INPUTDIR=/eos/uscms/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Zmumu_pT20
    OUTPUTDIR=$INPUTDIR/merged

    if [ ! -d $OUTPUTDIR ]; then
        echo "Creating output directory $OUTPUTDIR"
        mkdir $OUTPUTDIR
    fi

    for i in $(ls $INPUTDIR | grep root | awk -F'_' '{print $1}' | sort | uniq); do
        echo $i
        echo "hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root"
        hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root
    done
fi

if [ $doZee = True ]; then
    INPUTDIR=/eos/uscms/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Zee_pT20
    OUTPUTDIR=$INPUTDIR/merged

    if [ ! -d $OUTPUTDIR ]; then
        echo "Creating output directory $OUTPUTDIR"
        mkdir $OUTPUTDIR
    fi

    for i in $(ls $INPUTDIR | grep root | awk -F'_' '{print $1}' | sort | uniq); do
        echo $i
        echo "hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root"
        hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root
    done
fi

if [ $doWm = True ]; then
    INPUTDIR=/eos/uscms/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Wmunu
    OUTPUTDIR=$INPUTDIR/merged

    if [ ! -d $OUTPUTDIR ]; then
        echo "Creating output directory $OUTPUTDIR"
        mkdir $OUTPUTDIR
    fi

    for i in $(ls $INPUTDIR | grep "wm\|data" | grep root | awk -F'_' '{print $1}' | sort | uniq); do
        echo $i
        echo "hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root"
        hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root
    done
fi

if [ $doWe = True ]; then
    INPUTDIR=/eos/uscms/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Wenu
    OUTPUTDIR=$INPUTDIR/merged

    if [ ! -d $OUTPUTDIR ]; then
        echo "Creating output directory $OUTPUTDIR"
        mkdir $OUTPUTDIR
    fi

    for i in $(ls $INPUTDIR | grep "we\|data" | grep root | awk -F'_' '{print $1}' | sort | uniq); do
        echo $i
        echo "hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root"
        hadd $OUTPUTDIR/$i.root $INPUTDIR/$i*root
    done
fi
