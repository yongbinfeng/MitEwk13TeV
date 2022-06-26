#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

POSTFIX=
NTUPLEDIR_13=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/13TeV/probes
OUTPUTDIR_13=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/13TeV/results

# integrated luminosity for data
LUMI13=200.9


#
# Muon efficiencies
#############################################
## 13 TeV ###################################
#############################################
# trigger efficiency
root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuHLTEff_aMCxPythia/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuHLTEff_aMCxPythia/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)

root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuHLTEff_aMCxPythia/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuHLTEff_aMCxPythia/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuHLTEff/probes.root\"\)


# selection, etc
root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)

root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)

# standalone
root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)

root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)

#
# systematics
#
# # id+iso efficiency
# # # powheg x pythia
root -l -b -q plotEff.C+\(\"mu_sit.bins\",5,1,5,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
# # powheg x photos
root -l -b -q plotEff.C+\(\"mu_sit.bins\",6,1,6,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
# # # minlo 
# root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_minloxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_minlo/probes.root\"\)
# # # # Alternate BKG fit
root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,7,2,7,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff/probes.root\"\)
# # # tag pt cut
root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia_tagPt/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt/probes.root\"\)
## tag Max pt cut
#root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia_tagPt_Max/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt_Max/probes.root\"\)
#root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia_tagPt_Max${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPt_Max/probes.root\"\)
### tag Pt cut at 20
#root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPtMin/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia_tagPtMin/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPtMin/probes.root\"\)
#root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff_tagPtMin/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia_tagPtMin${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPtMin/probes.root\"\)
### tag Pt cut at 20-25GeV
#root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPtMin20_Max25/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuSITEff_aMCxPythia_tagPtMin20_Max25/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPtMin20_Max25/probes.root\"\)
#root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_13}/Zmm/Data/MuSITEff_tagPtMin20_Max25/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuSITEff_aMCxPythia_tagPtMin20_Max25${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuSITEff_tagPtMin20_Max25/probes.root\"\)



# # # standalone efficiency
#mc
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# # # # powheg x pythia
root -l -b -q plotEff.C+\(\"mu_sta.bins\",5,1,7,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# # # # # powheg x photos
root -l -b -q plotEff.C+\(\"mu_sta.bins\",6,1,7,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
# # # # # minlo
# root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,6,2,6,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_minloxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_minlo/probes.root\"\)
# # # bkg power law
root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,7,7,7,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff/probes.root\"\)
## tag pt cut
root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia_tagPt/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff_tagPt/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt/probes.root\"\)
## tag max pt cut
#root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia_tagPt_Max/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt_Max/probes.root\"\)
#root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia_tagPt_Max${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPt_Max/probes.root\"\)
### tag pt cut at 20
#root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPtMin/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia_tagPtMin/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPtMin/probes.root\"\)
#root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff_tagPtMin/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia_tagPtMin${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPtMin/probes.root\"\)
### tag pt cut at 20-25GeV
#root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPtMin20_Max25/probes.root\",\"${OUTPUTDIR_13}/Zmm/MC/MuStaEff_aMCxPythia_tagPtMin20_Max25/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPtMin20_Max25/probes.root\"\)
#root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_13}/Zmm/Data/MuStaEff_tagPtMin20_Max25/probes.root\",\"${OUTPUTDIR_13}/Zmm/Data/MuStaEff_aMCxPythia_tagPtMin20_Max25${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI13},\"${NTUPLEDIR_13}/Zmm/MC/MuStaEff_tagPtMin20_Max25/probes.root\"\)



# # ## Make a bunch of plots ##
TYPE=_aMCxPythia_tagPt
# # # ## TYPE
 # # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}NewTAble/pt/Negative\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}/Negative/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}/Negative/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mufineptbins.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/pt/Negative\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/pt/Positive\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/eta/Negative\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmHLT${TYPE}${POSTFIX}/eta/Positive\",\"$OUTPUTDIR_13/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_hlt.bins\"\)
  
 # # ## SEL
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSIT${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_sit.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSIT${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_sit.bins\"\)
 
 # # # # # # STA
 # root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSta${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.00,1.20,$LUMI13,\"mu_sta.bins\"\)
 
 # root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_13/Plots/ZmmSta${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_13/Zmm/MC/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_13/Zmm/Data/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI13,\"mu_sta.bins\"\)

# rm *.so *.d
