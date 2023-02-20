"""
plot and save the MET x and y corrections
for Z data and MC
"""

lumi = 200.84
fnames = {}

filedir = "/eos/home-y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/"
fnames["13TeV"] = {}
fnames["13TeV"]["ee"] = {}
fnames["13TeV"]["ee"]["Data"] = filedir + "Zee_pT20/merged/data.root"
fnames["13TeV"]["ee"]["MC"] = filedir + "Zee_pT20/merged/zee.root"
fnames["13TeV"]["mumu"] = {}
fnames["13TeV"]["mumu"]["Data"] = filedir + "Zmumu_pT20/merged/data.root"
fnames["13TeV"]["mumu"]["MC"] = filedir + "Zmumu_pT20/merged/zmm.root"

filedir = "/eos/home-y/yofeng/Ntuples_LowPU/5TeV/20230213/Selections/"
fnames["5TeV"] = {}
fnames["5TeV"]["ee"] = {}
fnames["5TeV"]["ee"]["Data"] = filedir + "Zee_pT20/merged/data.root"
fnames["5TeV"]["ee"]["MC"] = filedir + "Zee_pT20/merged/zee.root"
fnames["5TeV"]["mumu"] = {}
fnames["5TeV"]["mumu"]["Data"] = filedir + "Zmumu_pT20/merged/data.root"
fnames["5TeV"]["mumu"]["MC"] = filedir + "Zmumu_pT20/merged/zmm.root"

import ROOT
ROOT.gROOT.SetBatch(True)

for era in ["13TeV", "5TeV"]:
    for channel in ["mumu", "ee"]:
        fMC = ROOT.TFile(fnames[era][channel]["MC"])
        nEvtMC = fMC.Get("hGenWeights").Integral()
        tMC = fMC.Get("Events")

        fData = ROOT.TFile(fnames[era][channel]["Data"])
        tData = fData.Get("Events")

        hMET_x_MC = ROOT.TH1D("hMET_x_MC", "h0", 300, -50.0, 50)
        hMET_y_MC = ROOT.TH1D("hMET_y_MC", "h0", 300, -50.0, 50)

        hMET_x_Data = ROOT.TH1D("hMET_x_Data", "h0", 300, -50.0, 50)
        hMET_y_Data = ROOT.TH1D("hMET_y_Data", "h0", 300, -50.0, 50)

        tMC.Draw("met * TMath::Cos(metPhi) >> %s" % hMET_x_MC.GetName(), "(lep1.Pt() > 25.0 && fabs(lep1.Eta()) < 2.4 && lep2.Pt() > 25.0 && fabs(lep2.Eta())<2.4) * scale1fb * %f / %f" % (lumi, nEvtMC))
        tMC.Draw("met * TMath::Sin(metPhi) >> %s" % hMET_y_MC.GetName(), "(lep1.Pt() > 25.0 && fabs(lep1.Eta()) < 2.4 && lep2.Pt() > 25.0 && fabs(lep2.Eta())<2.4) * scale1fb * %f / %f" % (lumi, nEvtMC))

        tData.Draw("met * TMath::Cos(metPhi) >> %s" % hMET_x_Data.GetName(), "(lep1.Pt() > 25.0 && fabs(lep1.Eta()) < 2.4 && lep2.Pt() > 25.0 && fabs(lep2.Eta())<2.4)")
        tData.Draw("met * TMath::Sin(metPhi) >> %s" % hMET_y_Data.GetName(), "(lep1.Pt() > 25.0 && fabs(lep1.Eta()) < 2.4 && lep2.Pt() > 25.0 && fabs(lep2.Eta())<2.4)")
        
        print("data era %s, channel %s" % (era, channel))
        print("MC MET_x mean: %f" % hMET_x_MC.GetMean())
        print("MC MET_y mean: %f" % hMET_y_MC.GetMean())
        print("Data MET_x mean: %f" % hMET_x_Data.GetMean())
        print("Data MET_y mean: %f" % hMET_y_Data.GetMean())

        c1 = ROOT.TCanvas("c1_mc", "c1", 800, 600)
        c1.SetLogy()
        hMET_x_MC.SetLineColor(2)
        hMET_x_MC.SetMarkerColor(2)
        hMET_x_MC.Draw()
        hMET_y_MC.SetLineColor(4)
        hMET_y_MC.SetMarkerColor(4)
        hMET_y_MC.Draw("same")
        hMET_x_Data.SetLineColor(6)
        hMET_x_Data.SetMarkerColor(6)
        hMET_x_Data.Draw("same")
        hMET_y_Data.SetLineColor(8)
        hMET_y_Data.SetMarkerColor(8)
        hMET_y_Data.Draw("same")
        c1.Print("plots/met_xy_%s_%s.pdf" % (era, channel))
        
        # save output
        ofile = ROOT.TFile("data/met_xy_%s_%s.root" % (era, channel), "recreate")
        hMET_x_MC.SetDirectory(ofile)
        hMET_x_MC.Write()
        hMET_y_MC.SetDirectory(ofile)
        hMET_y_MC.Write()
        hMET_x_Data.SetDirectory(ofile)
        hMET_x_Data.Write()
        hMET_y_Data.SetDirectory(ofile)
        hMET_y_Data.Write()
        ofile.Write()
        ofile.Close()
        
        print("\n\n")