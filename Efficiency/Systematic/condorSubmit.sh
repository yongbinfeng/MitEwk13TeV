#!/bin/bash

exe=./runStep2Step3.sh
PU=`echo "${2}"`
ID=`echo "${3}"`
# NBINS=96 #96 bins for muon channels
# NBINS=24 #24 bins for STANDALONE MUON channels
NBINS=48 #for MuSITEff and EleGSFSelEff
NBINS=18 # for MuStaEff
#NTOYS=300 # should be default 1000
#EFFTYPE=MuSITEff
#EFFTYPE=MuStaEff
EFFTYPE=EleGSFSelEff

# declare -a EFFTYPES=("MuHLTEff" "MuSelEff" "MuStaEff") #for muons
declare -a CHARGES=("Combined") #for muons
# declare -a CHARGES=("Positive" "Negative") #for muons
#FOLDER=Zmm # or Zee
FOLDER=Zee # or Zee
 
VERS=

postfixs=(_POWxPythia${VERS} _aMCxPythia${VERS} _aMCxPythia${VERS})
postfixs_alt=(_POWxPhotos${VERS} _POWxPythia${VERS} _POWBKG${VERS})

# ## FSR alternative shape evaluation
#POSTFIX=_POWxPythia${VERS}
#POSTFIX_alt=_POWxPhotos${VERS}

# MC alternative shape eval
#POSTFIX=_aMCxPythia${VERS}
#POSTFIX_alt=_POWxPythia${VERS}

## BKG alternative shape eval
#POSTFIX=_aMCxPythia${VERS}
#POSTFIX_alt=_POWBKG${VERS}

for i in "${!postfixs[@]}"; do
    POSTFIX=${postfixs[i]}
    POSTFIX_alt=${postfixs_alt[i]}
    echo $POSTFIX
    echo $POSTFIX_alt
    echo "***"
    tmpname=tmp_${EFFTYPE}_${i}.sub
    echo $tmpname
    for CHARGE in "${CHARGES[@]}"
    do
        echo "submitting jobs to calculate efficiencies for ${EFFTYPE} ${CHARGE}, total bins ${NBINS} with ${NTOYS} toys" 
        echo "executable              = "${exe} > $tmpname
        echo "arguments               = \$(ClusterId) \$(ProcId)" ${NTOYS} ${FOLDER} ${EFFTYPE} ${CHARGE} ${POSTFIX} ${POSTFIX_alt}>> $tmpname
        echo "output                  = $PWD/output/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).out" >> $tmpname
        echo "error                   = $PWD/error/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).err"  >> $tmpname
        echo "log                     = $PWD/log/${EFFTYPE}.${CHARGE}.${NTOYS}.\$(ProcId).\$(ClusterId).log"    >> $tmpname
        echo "+JobFlavour = \"workday\" ">> $tmpname
        echo "queue ${NBINS}" >> $tmpname
        condor_submit $tmpname
    done
done
