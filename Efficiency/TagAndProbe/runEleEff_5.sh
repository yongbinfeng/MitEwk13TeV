#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

POSTFIX=
NTUPLEDIR_5=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/probes
OUTPUTDIR_5=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/results

# LUMI only used for plotting
LUMI5=298.0

#
# Electron efficiencies
#
# # # hlt efficiency
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/MC/EleHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/MC/EleHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zee/MC/EleHLTEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/Data/EleHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}/Negative\",\"png\",0,0,-1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/Data/EleHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleHLTEff_aMCxPythia${POSTFIX}/Positive\",\"png\",0,0,1,\"Supercluster\",\"supercluster\",0.0,1.15,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleHLTEff/probes.root\"\)

# # # # ####################### gsf+is+iso efficiency
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_5}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,2,2,2,\"${NTUPLEDIR_5}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff/probes.root\"\)

#
# gsf+is+iso systematics
#
### Tag Pt
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_5}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,2,2,2,\"${NTUPLEDIR_5}/Zee/Data/EleGSFSelEff_tagPt/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff_tagPt/probes.root\"\)
### tag max pt
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_5}/Zee/MC/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt_Max/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff_tagPt_Max/probes.root\"\)
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,2,2,2,\"${NTUPLEDIR_5}/Zee/Data/EleGSFSelEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleGSFSelEff_aMCxPythia${POSTFIX}_tagPt_Max/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff_tagPt_Max/probes.root\"\)
### bkg systematics
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",2,7,2,7,\"${NTUPLEDIR_5}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleGSFSelEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff/probes.root\"\)
### powheg x pythia
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",5,2,5,2,\"${NTUPLEDIR_5}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleGSFSelEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff/probes.root\"\)
### powheg x photos
root -l -b -q plotEff.C+\(\"ele_gsf.bins\",6,2,6,2,\"${NTUPLEDIR_5}/Zee/Data/EleGSFSelEff/probes.root\",\"${OUTPUTDIR_5}/Zee/Data/EleGSFSelEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Supercluster\",\"GSF\",0.5,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zee/MC/EleGSFSelEff/probes.root\"\)


TYPE=_aMCxPythia
## # # # #### Draw Plots
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_5/Plots/ZeeHLT${TYPE}${POSTFIX}/pt/Negative\",\"$OUTPUTDIR_5/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"PtBins\",0.01,1.20,0.7,1.1,$LUMI5,\"ele_hlt.bins\"\)
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_5/Plots/ZeeHLT${TYPE}${POSTFIX}/pt/Positive\",\"$OUTPUTDIR_5/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"PtBins\",0.01,1.20,0.7,1.1,$LUMI5,\"ele_hlt.bins\"\)
#
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_5/Plots/ZeeHLT${TYPE}${POSTFIX}/eta/Negative\",\"$OUTPUTDIR_5/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaBins\",0.00,1.20,$LUMI5,\"ele_hlt.bins\"\)
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_5/Plots/ZeeHLT${TYPE}${POSTFIX}/eta/Positive\",\"$OUTPUTDIR_5/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaBins\",0.00,1.20,$LUMI5,\"ele_hlt.bins\"\)
# 
##root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_5/Plots/ZeeHLT${TYPE}${POSTFIX}/etapt/Negative\",\"$OUTPUTDIR_5/Zee/MC/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI5,\"ele_hlt.bins\"\)
##root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_5/Plots/ZeeHLT${TYPE}${POSTFIX}/etapt/Positive\",\"$OUTPUTDIR_5/Zee/MC/EleHLTEff${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI5,\"ele_hlt.bins\"\)
# 
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_5/Plots/ZeeGSFSel${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_5/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.01,1.20,0.7,1.1,$LUMI5,\"ele_gsf.bins\"\)
# 
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_5/Plots/ZeeGSFSel${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_5/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI5,\"ele_gsf.bins\"\)
 
#root -l -b -q plotDataMC.C+\(\"$OUTPUTDIR_5/Plots/ZeeGSFSel${TYPE}${POSTFIX}/etapt/Combined\",\"$OUTPUTDIR_5/Zee/MC/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_5/Zee/Data/EleGSFSelEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaPtBins\",0.00,1.20,$LUMI5,\"ele_gsf.bins\"\)


# rm *.so *.d
