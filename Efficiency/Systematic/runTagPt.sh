#!/bin/bash
ClusterID=$1
binnum=$2
NTOYS=$3
FOLDER=$4
EFFTYPE=$5
CHARGE=$6
POSTFIX=$7
POSTFIX_alt=$8
POSTFIX_alt2=$9

WORKDIR="/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/MitEwk13TeV"
FILEDIR=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/results/

CMSSW_BASE="/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/"
TOP="/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/"
# echo ${ClusterID}
echo "running job w Cluster ID ${ClusterID}, binnum ${binnum}, NTOYS ${NTOYS}, FOLDER ${FOLDER}, EFFTYPE ${EFFTYPE}, Charge ${CHARGE}"

## FSR alternative shape evaluation
# POSTFIX=_POWxPythia_v1
# POSTFIX_alt=_POWxPhotos_v1

## MC alternative shape eval
# POSTFIX=_aMCxPythia_v1
# POSTFIX_alt=_minloxPythia_v1

## BKG alternative shape eval
# POSTFIX=_aMCxPythia_v1
# POSTFIX_alt=_POWBKG_v1
#
BINVAR=etapt # probably don't need to change
OUTPUTDIR=${FILEDIR}/TOYS/${EFFTYPE}${POSTFIX}${POSTFIX_alt}/${CHARGE}/
DIR1=${FILEDIR}/${FOLDER}/Data/${EFFTYPE}${POSTFIX}/${CHARGE}/plots/
DIR2=${FILEDIR}/${FOLDER}/Data/${EFFTYPE}${POSTFIX_alt}/${CHARGE}/plots/
DIR3=${FILEDIR}/${FOLDER}/Data/${EFFTYPE}${POSTFIX_alt2}/${CHARGE}/plots/
cd $CMSSW_BASE
eval `scramv1 runtime -sh`
cd $TOP
mkdir -p ${OUTPUTDIR}
# mkdir -p ${TOP}/${EFFTYPE}/${CHARGE}/Step2Output/${POSTFIX}v${POSTFIX_alt}
#TOP="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"
# root -l toyGenAndPull.C+
root -l -b << EOF
gSystem->Load("${WORKDIR}/Efficiency/Systematic/toyGenAndPull_TagPt_C.so")
toyGenAndPull_TagPt("${DIR2}","${DIR1}","${BINVAR}_${binnum}","${OUTPUTDIR}","var_${BINVAR}_${binnum}",${NTOYS}, "${DIR3}")
.q
EOF
