#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

POSTFIX=
NTUPLEDIR_13=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/13TeV/probes
OUTPUTDIR_13=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/13TeV/results

# integrated luminosity for data
LUMI13=200.9


#################################################################
# Electron efficiencies
#
# # # hlt efficiency
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/Data/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/Data/EleHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleHLTEff/probes.root\"\)

# # # # ####################### gsf+is+iso efficiency
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,2,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)

#
# systematics
#
# bkg systematic
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
### Tag Pt
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,2,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
### tag max pt
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_13}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt_Max/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt_Max/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,2,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt_Max/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_tagPt_Max/probes.root\"\)
### powheg x pythia
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",5,2,5,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)
### powheg x photos
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",6,2,6,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff/probes.root\"\)

# root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,1,2,2,\"${NTUPLEDIR_13}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_13}/Zee/Data/EleGSFSelEff_minloxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zee/MC/EleGSFSelEff_minlo/probes.root\"\)


TYPE=_aMCxPythia
# # # # #### Draw Plots
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/pt/Negative\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"PtBins\",0.00,1.20,0.7,1.1,$LUMI13,\"ele_hlt.bins\"\)
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/pt/Positive\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"PtBins\",0.00,1.20,0.67,1.1,$LUMI13,\"ele_hlt.bins\"\)
# 
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/eta/Negative\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/eta/Positive\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
#  
## root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/etapt/Negative\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
## root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_13/Plots/ZeeHLT${TYPE}${POSTFIX}/etapt/Positive\",\"$OUTPUTDIR_13/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI13,\"ele_hlt.bins\"\)
# 
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZeeGSFSel${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_13/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.00,1.20,0.7,1.1,$LUMI13,\"ele_gsf.bins\"\)
#
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZeeGSFSel${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_13/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"ele_gsf.bins\"\)
# 
##root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_13/Plots/ZeeGSFSel${TYPE}${POSTFIX}/etapt/Combined\",\"$OUTPUTDIR_13/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI13,\"ele_gsf.bins\"\)
#
## rm *.so *.d
