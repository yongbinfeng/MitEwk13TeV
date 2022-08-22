import ROOT
import sys
sys.path.append("/uscms/home/yfeng/nobackup/WpT/Cards/TestCode/WpTAnalysis")
from CMSPLOTS.myFunction import DrawHistos


ROOT.gROOT.SetBatch(True)


ifile = ROOT.TFile("/uscms/home/yfeng/nobackup/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Utils/L1MuonPrefiringParametriations.root")

hnames = [
    "L1prefiring_muonparam_0.0To0.2_20172018",
    "L1prefiring_muonparam_0.2To0.3_20172018",
    "L1prefiring_muonparam_0.3To0.55_20172018",
    "L1prefiring_muonparam_0.55To0.83_20172018",
    "L1prefiring_muonparam_0.83To1.24_20172018",
    "L1prefiring_muonparam_1.24To1.4_20172018",
    "L1prefiring_muonparam_1.6To1.8_20172018",
    "L1prefiring_muonparam_1.8To2.1_20172018",
    "L1prefiring_muonparam_2.1To2.25_20172018",
    "L1prefiring_muonparam_2.25To2.4_20172018",
]

colors = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 11
]

labels = [
    "0.0 < |#eta| < 0.2",
    "0.2 < |#eta| < 0.3",
    "0.3 < |#eta| < 0.55",
    "0.55 < |#eta| < 0.83",
    "0.83 < |#eta| < 1.24",
    "1.24 < |#eta| < 1.6",
    "1.6 < |#eta| < 1.8",
    "1.8 < |#eta| < 2.1",
    "2.1 < |#eta| < 2.25",
    "2.25 < |#eta| < 2.4",    
]

histos = []
for idx, hname in enumerate(hnames):
    histo = ifile.Get(hname).Clone(hname + "_Cloned")
    histo.SetLineColor(colors[idx])
    histo.SetMarkerColor(colors[idx])
    histos.append( histo ) 


DrawHistos(histos, labels, 0, 40, "#mu p_{T} [GeV]", 0, 0.011, "Prefiring Prob.", "L1prefiring_muon_20172018", donormalize=False, showratio=False, dology=False, noCMS=True, legendPos=[0.2, 0.55, 0.5, 0.9], noLumi=True, nMaxDigits=3)
