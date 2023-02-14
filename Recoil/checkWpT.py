"""
plot the W pt distribution from different samples,
to check if it is smooth, and if there are lots of zeros
"""

filedir = "/eos/home-y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/Wmunu/merged/"
fname0 = filedir + "wm0.root"
fname1 = filedir + "wm1.root"
fname2 = filedir + "wm2.root"

import ROOT

f0 = ROOT.TFile(fname0)
f1 = ROOT.TFile(fname1)
f2 = ROOT.TFile(fname2)

n0 = f0.Get("hGenWeights").Integral()
n1 = f1.Get("hGenWeights").Integral()
n2 = f2.Get("hGenWeights").Integral()

t0 = f0.Get("Events")
t1 = f1.Get("Events")
t2 = f2.Get("Events")

vnames = {
    "genVPt": "genV.Pt()",
    "genLepPt": "genLep.Pt()",
    "genNuPt": "genNu.Pt()",
}

for vname in ["genVPt", "genLepPt", "genNuPt"]:
    h0 = ROOT.TH1D("h0_" + vname, "h0", 150, 0, 150)
    h1 = ROOT.TH1D("h1_" + vname, "h1", 150, 0, 150)
    h2 = ROOT.TH1D("h2_" + vname, "h2", 150, 0, 150)
    
    t0.Draw(vnames[vname] + ">>h0_" + vname, "(lep.Pt() > 25.0 && fabs(lep.Eta()) < 2.4) * scale1fb / %f" % n0)
    t1.Draw(vnames[vname] + ">>h1_" + vname, "(lep.Pt() > 25.0 && fabs(lep.Eta()) < 2.4) * scale1fb / %f" % n1)
    t2.Draw(vnames[vname] + ">>h2_" + vname, "(lep.Pt() > 25.0 && fabs(lep.Eta()) < 2.4) * scale1fb / %f" % n2)
    
    ht = h0.Clone("htotal_" + vname)
    ht.Add(h1)
    ht.Add(h2)
    
    c1 = ROOT.TCanvas("c1_" + vname, "c1", 800, 600)
    c1.SetLogy()
    ht.SetLineColor(ROOT.kBlack)
    ht.Draw()
    h0.SetLineColor(ROOT.kRed)
    h0.Draw("same")
    h1.SetLineColor(ROOT.kBlue)
    h1.Draw("same")
    h2.SetLineColor(ROOT.kGreen)
    h2.Draw("same")
    
    c1.Print(vname + ".pdf")