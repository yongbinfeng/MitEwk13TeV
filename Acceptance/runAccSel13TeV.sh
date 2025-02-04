#! /bin/bash


SYS=""
OUTDIR=/uscms/home/yfeng/nobackup/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Acceptance/ResultsSEL
EFFDIR="/uscms_data/d3/yfeng/WpT/Data/lowpu_dataNew/13TeV/results"

S13="13TeV"
S5="5TeV"

#### Efficiency Systematics Files
EFFSYSDIR="/uscms_data/d3/yfeng/WpT/Data/lowpu_dataNew/13TeV/results/Systematics"
Sta=${EFFSYSDIR}/SysUnc_MuStaEff.root
SIT=${EFFSYSDIR}/SysUnc_MuSITEff.root
GSF=${EFFSYSDIR}/SysUnc_EleGSFSelEff.root

# #
# #    W- > mu nu
# # #

## Wm positive
root -l -q computeAccSelWm.C+\(\"w13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_wmp_13TeV\",\"ACC_wmp_13_0\",1,0,\"${SIT}\",\"${Sta}\",1\)

## Wm negative
root -l -q computeAccSelWm.C+\(\"w13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_wmm_13TeV\",\"ACC_wmm_13_0\",-1,0,\"${SIT}\",\"${Sta}\",1\)

## Wm (charge combined)
root -l -q computeAccSelWm.C+\(\"w13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_wm_13TeV\",\"ACC_wm_13_0\",0,0,\"${SIT}\",\"${Sta}\",1\)

####
####     Z -> mu mu
root -l -q computeAccSelZmm.C+\(\"z13.conf\",\"${EFFDIR}/Zmm/\",\"${OUTDIR}/SEL_zmm_13TeV${SYS}\",\"ACC_zmm_13\",0,\"${SIT}\",\"${Sta}\",1\)


# #
# # ##    W -> e nu
# #

## We positive
root -l -q computeAccSelWe.C+\(\"w13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_wep_13TeV${SYS}\",\"ACC_wep_13_0\",1,0,1,0,\"${GSF}\",1\)

## We negative
root -l -q computeAccSelWe.C+\(\"w13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_wem_13TeV${SYS}\",\"ACC_wem_13_0\",-1,0,1,0,\"${GSF}\",1\)

## We (charge combined)
root -l -q computeAccSelWe.C+\(\"w13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_we_13TeV${SYS}\",\"ACC_we_13_0\",0,0,1,0,\"${GSF}\",1\)



# # #
# # # Z->ee
# # #
root -l -q computeAccSelZee.C+\(\"z13.conf\",\"${EFFDIR}/Zee/\",\"${OUTDIR}/SEL_zee_13TeV${SYS}\",\"ACC_zee_13\",\"${S13}\",0,1,0,\"${GSF}\",1\)


# rm *.so *.d
