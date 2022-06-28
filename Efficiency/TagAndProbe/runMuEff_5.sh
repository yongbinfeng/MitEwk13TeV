#--------------------------------------------------------------
#  Example script to run efficiency calculations
#==============================================================

POSTFIX=
NTUPLEDIR_5=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/probes
OUTPUTDIR_5=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/results
LUMI5=298.0


#
# Muon efficiencies
#

### --------------------------- HLT Eff MC
root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuHLTEff_aMCxPythia/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuHLTEff_aMCxPythia/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\"\)

### ---------------------------  HLT EFF Data
root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuHLTEff_aMCxPythia/Positive\",\"png\",0,0,1,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_hlt.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/Data/MuHLTEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuHLTEff_aMCxPythia/Negative\",\"png\",0,0,-1,\"Muon\",\"trigger\",0.7,1.03,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuHLTEff/probes.root\"\)

### --------------------------- SIT Eff MC
root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuSITEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)

### --------------------------- SIT Eff Data
root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSITEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)

### --------------------------- Standalone Eff MC
root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuStaEff_aMCxPythia/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)

### --------------------------- Standalone Eff Data
# change muon standalone bkg function to exponential function as well, orignal was 6 quadratic function: 2,6,2,6
root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuStaEff_aMCxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)


### ---------------------------- Selection ID Iso, alternative models
# # powheg x pythia: MC systematic
root -l -b -q plotEff.C+\(\"mu_sit.bins\",5,1,5,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSITEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)
# # powheg x photos: FSR systematic
root -l -b -q plotEff.C+\(\"mu_sit.bins\",6,1,6,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSITEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)
# # Alternate BKG fit: BKG systematic
root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,7,2,7,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSITEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff/probes.root\"\)
## tag pt cut
root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff_tagPt/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuSITEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff_tagPt/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff_tagPt/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSITEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff_tagPt/probes.root\"\)
## tag Max pt cut
root -l -b -q plotEff.C+\(\"mu_sit.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuSITEff_aMCxPythia_tagPt_Max/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff_tagPt_Max/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_sit.bins\",2,1,2,1,\"${NTUPLEDIR_5}/Zmm/Data/MuSITEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuSITEff_aMCxPythia_tagPt_Max${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"ID+Iso+Trk\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuSITEff_tagPt_Max/probes.root\"\)



### ---------------------------- Standalone, alternative models
# # powheg x pythia
root -l -b -q plotEff.C+\(\"mu_sta.bins\",5,1,7,1,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuStaEff_POWxPythia${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)
# # powheg x photos
root -l -b -q plotEff.C+\(\"mu_sta.bins\",6,1,7,1,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuStaEff_POWxPhotos${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)
# # bkg power law
root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,7,7,7,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuStaEff_POWBKG${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff/probes.root\"\)
## tag pt cut
root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_tagPt/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuStaEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_tagPt/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff_tagPt/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuStaEff_aMCxPythia_tagPt${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_tagPt/probes.root\"\)
### tag pt cut max
root -l -b -q plotEff.C+\(\"mu_sta.bins\",0,0,0,0,\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_5}/Zmm/MC/MuStaEff_aMCxPythia_tagPt_Max/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_tagPt_Max/probes.root\"\)
root -l -b -q plotEff.C+\(\"mu_sta.bins\",2,1,7,1,\"${NTUPLEDIR_5}/Zmm/Data/MuStaEff_tagPt_Max/probes.root\",\"${OUTPUTDIR_5}/Zmm/Data/MuStaEff_aMCxPythia_tagPt_Max${POSTFIX}/Combined\",\"png\",0,0,0,\"Muon\",\"stand-alone\",0.7,1.02,${LUMI5},\"${NTUPLEDIR_5}/Zmm/MC/MuStaEff_tagPt_Max/probes.root\"\)


#
# # ## Make a bunch of plots ##
#
TYPE=_aMCxPythia
## # # ## HLT
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_5/Plots/ZmmHLT${TYPE}${POSTFIX}/pt/Negative\",\"$OUTPUTDIR_5/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"PtBins\",0.01,1.20,0.9,1.1,$LUMI5,\"mu_hlt.bins\"\)
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_5/Plots/ZmmHLT${TYPE}${POSTFIX}/pt/Positive\",\"$OUTPUTDIR_5/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"PtBins\",0.01,1.20,0.9,1.1,$LUMI5,\"mu_hlt.bins\"\)
#
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_5/Plots/ZmmHLT${TYPE}${POSTFIX}/eta/Negative\",\"$OUTPUTDIR_5/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Negative/eff.root\",\"EtaBins\",0.00,1.20,$LUMI5,\"mu_hlt.bins\"\)
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_5/Plots/ZmmHLT${TYPE}${POSTFIX}/eta/Positive\",\"$OUTPUTDIR_5/Zmm/MC/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuHLTEff${TYPE}${POSTFIX}/Positive/eff.root\",\"EtaBins\",0.00,1.20,$LUMI5,\"mu_hlt.bins\"\)
#  
## # # ## SEL
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_5/Plots/ZmmSIT${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_5/Zmm/MC/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.01,1.20,0.96,1.04,$LUMI5,\"mu_sit.bins\"\)
#
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_5/Plots/ZmmSIT${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_5/Zmm/MC/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuSITEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI5,\"mu_sit.bins\"\)
# 
## # # # STA
root -l -b -q plotDataMC_singlepTbins.C+\(\"$OUTPUTDIR_5/Plots/ZmmSta${TYPE}${POSTFIX}/pt/Combined\",\"$OUTPUTDIR_5/Zmm/MC/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"PtBins\",0.01,1.20,0.96,1.04,$LUMI5,\"mu_sta.bins\"\)
# 
#root -l -b -q plotDataMC_singleEtabins.C+\(\"$OUTPUTDIR_5/Plots/ZmmSta${TYPE}${POSTFIX}/eta/Combined\",\"$OUTPUTDIR_5/Zmm/MC/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"$OUTPUTDIR_5/Zmm/Data/MuStaEff${TYPE}${POSTFIX}/Combined/eff.root\",\"EtaBins\",0.00,1.20,$LUMI5,\"mu_sta.bins\"\)

# rm *.so *.d
