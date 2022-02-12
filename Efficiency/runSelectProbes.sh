#! /bin/bash

INPUTDIR=/eos/home-y/yofeng/LowPU/Selection
OUTPUTDIR=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNewMeasured/13TeV/probes
#
# Select probes for muon efficiencies
#

# ## Zmm aMC@NLO with weights for the pythia and powheg w/ photos reweightings
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuHLTEff\",0,0,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuSITEff\",8,0,1\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff\",4,0,1\)

# ## Tag PT cut
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuHLTEff_tagPt\",0,0,1,30.0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuSITEff_tagPt\",8,0,1,30.0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff_tagPt\",4,0,1,30.0\)

### Zmm Minlo sample for one of the uncertainties
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_minlo/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuHLTEff_minlo\",0,0,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_minlo/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuSITEff_minlo\",8,0,1\)
#root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu_minlo/ntuples_0_1/zmm_select.raw.root\",\"${OUTPUTDIR}/Zmm/MC/MuStaEff_minlo\",4,0,1\)

# Zmm Data
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuHLTEff\",0,0,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuSITEff\",8,0,0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuStaEff\",4,0,0\)

## tag pt cut
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuHLTEff_tagPt\",0,0,0,30.0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuSITEff_tagPt\",8,0,0,30.0\)
root -l -n -q selectProbesMuEff.C+\(\"${INPUTDIR}/Zmumu/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zmm/Data/MuStaEff_tagPt\",4,0,0,30.0\)

# These are the ones for low PU 13 TeV
# mc
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleHLTEff\",0,0,1,0\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleGSFSelEff\",4,0,1,0\)

root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleHLTEff_tagPt\",0,0,1,0,30\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleGSFSelEff_tagPt\",4,0,1,0,30\)

#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_minlo/ntuples_0_1/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleHLTEff_minlo\",0,0,1,0\)
#root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee_minlo/ntuples_0_1/zee_select.root\",\"${OUTPUTDIR}/Zee/MC/EleGSFSelEff_minlo\",4,0,1,0\)
# # # # data
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleHLTEff\",0,0,0,0\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleGSFSelEff\",4,0,0,0\)

root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleHLTEff_tagPt\",0,0,0,0,30\)
root -l -n -q selectProbesEleEff.C+\(\"${INPUTDIR}/Zee/ntuples_0_1/data_select.root\",\"${OUTPUTDIR}/Zee/Data/EleGSFSelEff_tagPt\",4,0,0,0,30\)


# rm *.so *.d
