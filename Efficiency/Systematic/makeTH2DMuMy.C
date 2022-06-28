#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMinuit.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../../Utils/MitStyleRemix.hh" // style settings for drawing

TString mainDir = "/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/results/TOYS/";
TString sigDirFSR = "_POWxPythia_POWxPhotos/";
//TString sigDirMC = "_aMCxPythia_minloxPythia/";
TString sigDirMC = "_aMCxPythia_POWxPythia/";
TString bkgDir = "_aMCxPythia_POWBKG/"; // should be exp vs Powerlaw
TString tagPtDir = "_aMCxPythia_aMCxPythia_tagPt/";

TString effDirD = "/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/results/Zmm/Data/";
TString effDirM = "/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/results/Zmm/MC/";
TString fsr = "_POWxPhotos/";
TString mc = "_POWxPythia/";
TString bkg = "_POWBKG/";
TString amc = "_aMCxPythia/";
TString tagpt = "_aMCxPythia_tagPt/";

TString outDir = "/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/lowpu_dataNew/5TeV/results/Systematics";
TString subf = "";

const vector<TString> charges{ "Combined", "Combined" };

void makeHTML(TString suffix) {
    ofstream htmlfile;
    //char htmlfname[500];
    //sprintf(htmlfname, "plots/plots.html", outDir.Data());
    TString htmlfname = outDir + "/plots_" + suffix + "/plots.html";
    htmlfile.open(htmlfname.Data());
    htmlfile << "<!DOCTYPE html" << endl;
    htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
    htmlfile << "<html>" << endl;

    htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
    htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;

    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_amcPosDat.png\"><img src=\"effetapt_amcPosDat.png\" alt=\"effetapt_amcPosDat.png\" width=\"100%\"><figcaption>effetapt_Data</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_amcPosMC.png\"><img src=\"effetapt_amcPosMC.png\" alt=\"effetapt_amcPosMC.png\" width=\"100%\"><figcaption>effetapt_MC</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_amcCorr.png\"><img src=\"effetapt_amcCorr.png\" alt=\"effetapt_amcCorr.png\" width=\"100%\"><figcaption>effetapt_MC</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;

    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_mcPosDat.png\"><img src=\"effetapt_mcPosDat.png\" alt=\"effetapt_mcPosDat.png\" width=\"100%\"><figcaption>effetapt_mc_Data</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_mcPosDatRatio.png\"><img src=\"effetapt_mcPosDatRatio.png\" alt=\"effetapt_mcPosDatRatio.png\" width=\"100%\"><figcaption>effetapt_mc_DataRatio</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_mcPosSys.png\"><img src=\"effetapt_mcPosSys.png\" alt=\"effetapt_mcPosSys.png\" width=\"100%\"><figcaption>effetapt_mc_systematic</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;

    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_fsrPosDat.png\"><img src=\"effetapt_fsrPosDat.png\" alt=\"effetapt_fsrPosDat.png\" width=\"100%\"><figcaption>effetapt_fsr_Data</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_fsrPosDatRatio.png\"><img src=\"effetapt_fsrPosDatRatio.png\" alt=\"effetapt_fsrPosDatRatio.png\" width=\"100%\"><figcaption>effetapt_fsr_DataRatio</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_fsrPosSys.png\"><img src=\"effetapt_fsrPosSys.png\" alt=\"effetapt_fsrPosSys.png\" width=\"100%\"><figcaption>effetapt_fsr_systematic</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;

    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_bkgPosDat.png\"><img src=\"effetapt_bkgPosDat.png\" alt=\"effetapt_bkgPosDat.png\" width=\"100%\"><figcaption>effetapt_bkg_Data</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_bkgPosDatRatio.png\"><img src=\"effetapt_bkgPosDatRatio.png\" alt=\"effetapt_bkgPosDatRatio.png\" width=\"100%\"><figcaption>effetapt_bkg_DataRatio</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_bkgPosSys.png\"><img src=\"effetapt_bkgPosSys.png\" alt=\"effetapt_bkgPosSys.png\" width=\"100%\"><figcaption>effetapt_bkg_systematic</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;

    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_tagPosDat.png\"><img src=\"effetapt_tagPosDat.png\" alt=\"effetapt_tagPosDat.png\" width=\"100%\"><figcaption>effetapt_tagPt_Data</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_tagPosMC.png\"><img src=\"effetapt_tagPosMC.png\" alt=\"effetapt_tagPosMC.png\" width=\"100%\"><figcaption>effetapt_tagPt_MC</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_tagPosDatRatio.png\"><img src=\"effetapt_tagPosDatRatio.png\" alt=\"effetapt_tagPosDatRatio.png\" width=\"100%\"><figcaption>effetapt_tagPt_DatRatio</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_tagPosMCRatio.png\"><img src=\"effetapt_tagPosMCRatio.png\" alt=\"effetapt_tagPosMCRatio.png\" width=\"100%\"><figcaption>effetapt_tagPt_MCRatio</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;

    htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_tagPosDatMCRatio.png\"><img src=\"effetapt_tagPosDatMCRatio.png\" alt=\"effetapt_tagPosDatMCRatio.png\" width=\"100%\"><figcaption>effetapt_tagPt_DatMCRatio</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_tagPosToyDatRatio.png\"><img src=\"effetapt_tagPosToyDatRatio.png\" alt=\"effetapt_tagPosToyDatRatio.png\" width=\"100%\"><figcaption>effetapt_tagPt_ToyRatio</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"effetapt_tagPosSys.png\"><img src=\"effetapt_tagPosSys.png\" alt=\"effetapt_tagPosSys.png\" width=\"100%\"><figcaption>effetapt_tagPt_systematic</figcaption></a></td>" << endl;
    htmlfile << "<td width=\"25%\"></td>" << endl;
    htmlfile << "</tr>" << endl;

}

void makeTH2DMuMy(TString eType = "MuStaEff")
{
    float ptrange_sit[5] = {25.0, 30, 35, 40, 50};
    float etarange_sit[13] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, 0., 0.3, 0.9, 1.2, 1.6, 2.1, 2.4};
    float ptrange_sta[2] = {25.0, 50.0};
    float etarange_sta[19] = {-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, -0.15, 0., 0.15, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4};
    int NBpt, NBeta;
    float *ptrange = NULL;
    float *etarange = NULL;
    if (eType == "MuStaEff") {
        NBpt = 1;
        NBeta = 18;
        ptrange = ptrange_sta;
        etarange= etarange_sta;
    } else {
        NBpt = 4;
        NBeta = 12;
        ptrange = ptrange_sit;
        etarange = etarange_sit;
    }
    bool doAbs = false;
    vector<double> vMCNeg;
    vector<double> vMCPos;
    vector<double> vFSRNeg;
    vector<double> vFSRPos;
    vector<double> vBkgNeg;
    vector<double> vBkgPos;
    vector<double> vTagNeg;
    vector<double> vTagPos;

    double value = 0; // read in value is the % diff from the central value
    char infilename[200];
    // Data
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), amc.Data(), subf.Data(), charges[0].Data());
    TFile* amcPd = new TFile(infilename);
    TH2D* amcPosDat = (TH2D*)amcPd->Get("hEffEtaPt");
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), amc.Data(), subf.Data(), charges[1].Data());
    TFile* amcNd = new TFile(infilename);
    TH2D* amcNegDat = (TH2D*)amcNd->Get("hEffEtaPt");
    // MC
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirM.Data(), eType.Data(), amc.Data(), subf.Data(), charges[0].Data());
    TFile* amcPm = new TFile(infilename);
    TH2D* amcPosMC = (TH2D*)amcPm->Get("hEffEtaPt");
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirM.Data(), eType.Data(), amc.Data(), subf.Data(), charges[1].Data());
    TFile* amcNm = new TFile(infilename);
    TH2D* amcNegMC = (TH2D*)amcNm->Get("hEffEtaPt");

    // fsr systematic
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), fsr.Data(), subf.Data(), charges[0].Data());
    TFile* fsrPd = new TFile(infilename);
    TH2D* fsrPosDat = (TH2D*)fsrPd->Get("hEffEtaPt");
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), fsr.Data(), subf.Data(), charges[1].Data());
    TFile* fsrNd = new TFile(infilename);
    TH2D* fsrNegDat = (TH2D*)fsrNd->Get("hEffEtaPt");

    //  MC systematic
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), mc.Data(), subf.Data(), charges[0].Data());
    TFile* mcPd = new TFile(infilename);
    TH2D* mcPosDat = (TH2D*)mcPd->Get("hEffEtaPt");
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), mc.Data(), subf.Data(), charges[1].Data());
    TFile* mcNd = new TFile(infilename);
    TH2D* mcNegDat = (TH2D*)mcNd->Get("hEffEtaPt");

    // bkg systematic
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), bkg.Data(), subf.Data(), charges[0].Data());
    TFile* bkgPd = new TFile(infilename);
    TH2D* bkgPosDat = (TH2D*)bkgPd->Get("hEffEtaPt");
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), bkg.Data(), subf.Data(), charges[1].Data());
    TFile* bkgNd = new TFile(infilename);
    TH2D* bkgNegDat = (TH2D*)bkgNd->Get("hEffEtaPt");

    // tag systematic
    // Data
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), tagpt.Data(), subf.Data(), charges[0].Data());
    TFile* tagPd = new TFile(infilename);
    TH2D* tagPosDat = (TH2D*)tagPd->Get("hEffEtaPt");
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirD.Data(), eType.Data(), tagpt.Data(), subf.Data(), charges[1].Data());
    TFile* tagNd = new TFile(infilename);
    TH2D* tagNegDat = (TH2D*)tagNd->Get("hEffEtaPt");
    // MC
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirM.Data(), eType.Data(), tagpt.Data(), subf.Data(), charges[0].Data());
    TFile* tagPm = new TFile(infilename);
    std::cout << "infilename name " << infilename << std::endl;
    TH2D* tagPosMC = (TH2D*)tagPm->Get("hEffEtaPt");
    std::cout <<  "******** " << tagPosMC->GetBinContent(1,1) << std::endl;
    sprintf(infilename, "%s%s%s%s%s/eff.root", effDirM.Data(), eType.Data(), tagpt.Data(), subf.Data(), charges[1].Data());
    TFile* tagNm = new TFile(infilename);
    TH2D* tagNegMC = (TH2D*)tagNm->Get("hEffEtaPt");

    for (int i = 0; i < NBeta * NBpt; ++i) {

        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), sigDirFSR.Data(), subf.Data(), charges[0].Data(), i);
        ifstream infile1(infilename);
        assert(infile1);
        infile1 >> value;
        value = (doAbs ? fabs(value) : value);
        vFSRNeg.push_back(value);
        value = 0;
        //std::cout << infilename << " "<< vFSRNeg.back() << std::endl;
        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), sigDirFSR.Data(), subf.Data(), charges[1].Data(), i);
        ifstream infile2(infilename);
        assert(infile2);
        infile2 >> value;
        value = (doAbs ? fabs(value) : value);
        vFSRPos.push_back(value);
        value = 0;
        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), sigDirMC.Data(), subf.Data(), charges[0].Data(), i);
        ifstream infile3(infilename);
        assert(infile3);
        infile3 >> value;
        value = (doAbs ? fabs(value) : value);
        vMCNeg.push_back(value);
        value = 0;
        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), sigDirMC.Data(), subf.Data(), charges[1].Data(), i);
        ifstream infile4(infilename);
        assert(infile4);
        infile4 >> value;
        value = (doAbs ? fabs(value) : value);
        vMCPos.push_back(value);
        value = 0;
        //std::cout << infilename << vMCPos.back() << std::endl;
        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), bkgDir.Data(), subf.Data(), charges[0].Data(), i);
        ifstream infile5(infilename);
        assert(infile5);
        infile5 >> value;
        value = (doAbs ? fabs(value) : value);
        vBkgNeg.push_back(value);
        value = 0;
        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), bkgDir.Data(), subf.Data(), charges[1].Data(), i);
        ifstream infile6(infilename);
        infile6 >> value;
        value = (doAbs ? fabs(value) : value);
        vBkgPos.push_back(value);
        value = 0;

        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), tagPtDir.Data(), subf.Data(), charges[0].Data(), i);
        ifstream infile7(infilename);
        assert(infile7);
        infile7 >> value;
        value = (doAbs ? fabs(value) : value);
        vTagNeg.push_back(value);
        value = 0;
        std::cout << infilename << vTagNeg.back() << std::endl;
        sprintf(infilename, "%s%s%s%s%s/var_etapt_%d.txt", mainDir.Data(), eType.Data(), tagPtDir.Data(), subf.Data(), charges[1].Data(), i);
        ifstream infile8(infilename);
        infile8 >> value;
        value = (doAbs ? fabs(value) : value);
        vTagPos.push_back(value);
        value = 0;
    }

    char histname[100];
    sprintf(histname, "hMCNeg");
    TH2D* hMCNeg = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hMCPos");
    TH2D* hMCPos = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hFSRNeg");
    TH2D* hFSRNeg = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hFSRPos");
    TH2D* hFSRPos = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hBkgNeg");
    TH2D* hBkgNeg = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hBkgPos");
    TH2D* hBkgPos = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    
    // the efficiencies dervied using a different tag-pt cut
    sprintf(histname, "hEffTagNeg");
    TH2D* hEffTagNeg = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hEffTagPos");
    TH2D* hEffTagPos = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hTagNeg");
    TH2D* hTagNeg = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);
    sprintf(histname, "hTagPos");
    TH2D* hTagPos = new TH2D(histname, histname, NBeta, etarange, NBpt, ptrange);

    for (int ipt = 0; ipt < NBpt; ipt++) {

        for (int ieta = 0; ieta < NBeta; ieta++) {

            // std::cout << ipt*NBeta+ieta << "  " << 1.+vFSRNeg[ipt*NBeta+ieta] << " " << 1.+vMCNeg[ipt*NBeta+ieta] << std::endl;
            double val_pos = amcPosDat->GetBinContent(ieta + 1, ipt + 1);
            double val_neg = amcNegDat->GetBinContent(ieta + 1, ipt + 1);
            // vector saves the absolute eff_default - eff_alternative
            // so sys_alternative = (eff_default - diff) / eff_default = 1.0 - diff / eff_default
            hFSRNeg->SetBinContent(ieta + 1, ipt + 1, 1. - vFSRNeg[ipt * NBeta + ieta] / val_neg);
            hFSRPos->SetBinContent(ieta + 1, ipt + 1, 1. - vFSRPos[ipt * NBeta + ieta] / val_pos);
            hMCNeg->SetBinContent(ieta + 1, ipt + 1,  1. - vMCNeg[ipt * NBeta + ieta] / val_neg);
            hMCPos->SetBinContent(ieta + 1, ipt + 1,  1. - vMCPos[ipt * NBeta + ieta] / val_pos);
            hBkgNeg->SetBinContent(ieta + 1, ipt + 1, 1. - vBkgNeg[ipt * NBeta + ieta] / val_neg);
            hBkgPos->SetBinContent(ieta + 1, ipt + 1, 1. - vBkgPos[ipt * NBeta + ieta] / val_pos);

            //double eff_fsr_pos = fsrPosDat->GetBinContent(ieta + 1, ipt + 1) / amcPosDat->GetBinContent(ieta + 1, ipt + 1);
            //hFSRPos->SetBinContent(ieta + 1, ipt + 1, eff_fsr_pos);
            //double eff_fsr_neg = fsrNegDat->GetBinContent(ieta + 1, ipt + 1) / amcNegDat->GetBinContent(ieta + 1, ipt + 1);
            //hFSRNeg->SetBinContent(ieta + 1, ipt + 1, eff_fsr_neg);

            //double eff_mc_pos = mcPosDat->GetBinContent(ieta + 1, ipt + 1) / amcPosDat->GetBinContent(ieta + 1, ipt + 1);
            //hMCPos->SetBinContent(ieta + 1, ipt + 1, eff_mc_pos);
            //double eff_mc_neg = mcNegDat->GetBinContent(ieta +1, ipt + 1) / amcNegDat->GetBinContent(ieta + 1, ipt + 1);
            //hMCNeg->SetBinContent(ieta + 1, ipt + 1, eff_mc_neg);

            //double eff_bkg_pos = bkgPosDat->GetBinContent(ieta + 1, ipt + 1) / amcPosDat->GetBinContent(ieta + 1, ipt + 1);
            //hBkgPos->SetBinContent(ieta + 1, ipt + 1, eff_bkg_pos);
            //double eff_bkg_neg = bkgNegDat->GetBinContent(ieta + 1, ipt + 1) / amcNegDat->GetBinContent(ieta + 1, ipt + 1);
            //hBkgNeg->SetBinContent(ieta + 1, ipt + 1, eff_bkg_neg);
            
            // tagPt systematic
            // vTag saves the absolute eff_default - eff_tagPt
            // so eff_tagPt = eff_default - diff
            hEffTagNeg->SetBinContent(ieta + 1, ipt + 1, amcNegDat->GetBinContent(ieta + 1, ipt + 1)  - vTagNeg[ipt * NBeta + ieta]);
            double num = 1.0 - vTagNeg[ipt * NBeta + ieta] / amcNegDat->GetBinContent(ieta + 1, ipt + 1);
            double dnm = tagNegMC->GetBinContent(ieta + 1, ipt + 1) / amcPosMC->GetBinContent(ieta + 1, ipt + 1);
            // if(num < 0.99 )num = 0.995;
            hTagNeg->SetBinContent(ieta + 1, ipt + 1, num / dnm);

            hEffTagPos->SetBinContent(ieta + 1, ipt + 1,  amcPosDat->GetBinContent(ieta + 1, ipt + 1) - vTagPos[ipt * NBeta + ieta]);
            num = 1.0 - vTagPos[ipt * NBeta + ieta] / amcPosDat->GetBinContent(ieta + 1, ipt + 1);
            dnm = tagPosMC->GetBinContent(ieta + 1, ipt + 1) / amcPosMC->GetBinContent(ieta + 1, ipt + 1);
            // if(num < 0.99 )num = 0.995;
            hTagPos->SetBinContent(ieta + 1, ipt + 1, num / dnm);
        }
    }

    // plus and minus are probably the same, so only compare one difference
    gSystem->mkdir((outDir+"/plots_"+eType).Data(), kTRUE);

    gStyle->SetPaintTextFormat(".3f");
    TCanvas* c = MakeCanvas("c", "c", 800, 600);
    //gStyle->SetPalette(1);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);

    CPlot::sOutDir = "";

    bool doLogy = 0;
    CPlot plotEffEtaPt_amcPosDat((outDir+"/plots_"+eType+"/effetapt_amcPosDat").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_amcPosDat.SetLogy(doLogy);
    plotEffEtaPt_amcPosDat.AddHist2D(amcPosDat, "COLZ,text");
    plotEffEtaPt_amcPosDat.Draw(c, 1, "png");

    CPlot plotEffEtaPt_amcPosMC((outDir+"/plots_"+eType+"/effetapt_amcPosMC").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_amcPosMC.SetLogy(doLogy);
    plotEffEtaPt_amcPosMC.AddHist2D(amcPosMC, "COLZ,text");
    plotEffEtaPt_amcPosMC.Draw(c, 1, "png");

    // fsr
    // efficiency using alternative FSR fit on data
    CPlot plotEffEtaPt_fsrPosDat((outDir+"/plots_"+eType+"/effetapt_fsrPosDat").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_fsrPosDat.SetLogy(doLogy);
    plotEffEtaPt_fsrPosDat.AddHist2D(fsrPosDat, "COLZ,text");
    plotEffEtaPt_fsrPosDat.Draw(c, 1, "png");
    // ratio betwee data fit using alternative model and default model
    fsrPosDat->Divide(amcPosDat);
    CPlot plotEffEtaPt_fsrPosDatRatio((outDir+"/plots_"+eType+"/effetapt_fsrPosDatRatio").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_fsrPosDatRatio.SetLogy(doLogy);
    plotEffEtaPt_fsrPosDatRatio.AddHist2D(fsrPosDat, "COLZ,text");
    plotEffEtaPt_fsrPosDatRatio.Draw(c, 1, "png");
    // fsr systematics
    CPlot plotEffEtaPt_fsrPosSys((outDir+"/plots_"+eType+"/effetapt_fsrPosSys").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_fsrPosSys.SetLogy(doLogy);
    plotEffEtaPt_fsrPosSys.AddHist2D(hFSRPos, "COLZ,text");
    plotEffEtaPt_fsrPosSys.Draw(c, 1, "png");

    // MC
    // efficiency using alternative MC fit on data
    CPlot plotEffEtaPt_mcPosDat((outDir+"/plots_"+eType+"/effetapt_mcPosDat").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_mcPosDat.SetLogy(doLogy);
    plotEffEtaPt_mcPosDat.AddHist2D(mcPosDat, "COLZ,text");
    plotEffEtaPt_mcPosDat.Draw(c, 1, "png");
    // ratio between data fit using alternative model and default model
    mcPosDat->Divide(amcPosDat);
    CPlot plotEffEtaPt_mcPosDatRatio((outDir+"/plots_"+eType+"/effetapt_mcPosDatRatio").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_mcPosDatRatio.SetLogy(doLogy);
    plotEffEtaPt_mcPosDatRatio.AddHist2D(mcPosDat, "COLZ,text");
    plotEffEtaPt_mcPosDatRatio.Draw(c, 1, "png");
    // MC systematic
    CPlot plotEffEtaPt_mcPosSys((outDir+"/plots_"+eType+"/effetapt_mcPosSys").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_mcPosSys.SetLogy(doLogy);
    plotEffEtaPt_mcPosSys.AddHist2D(hMCPos, "COLZ,text");
    plotEffEtaPt_mcPosSys.Draw(c, 1, "png");

    // bkg
    // efficiency using alternative bkg model to fit on data
    CPlot plotEffEtaPt_bkgPosDat((outDir+"/plots_"+eType+"/effetapt_bkgPosDat").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_bkgPosDat.SetLogy(doLogy);
    plotEffEtaPt_bkgPosDat.AddHist2D(bkgPosDat, "COLZ,text");
    plotEffEtaPt_bkgPosDat.Draw(c, 1, "png");
    // ratio between alternative fit and default model
    bkgPosDat->Divide(amcPosDat);
    CPlot plotEffEtaPt_bkgPosDatRatio((outDir+"/plots_"+eType+"/effetapt_bkgPosDatRatio").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_bkgPosDatRatio.SetLogy(doLogy);
    plotEffEtaPt_bkgPosDatRatio.AddHist2D(bkgPosDat, "COLZ,text");
    plotEffEtaPt_bkgPosDatRatio.Draw(c, 1, "png");
    // bkg systematics
    CPlot plotEffEtaPt_bkgPosSys((outDir+"/plots_"+eType+"/effetapt_bkgPosSys").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_bkgPosSys.SetLogy(doLogy);
    plotEffEtaPt_bkgPosSys.AddHist2D(hBkgPos, "COLZ,text");
    plotEffEtaPt_bkgPosSys.Draw(c, 1, "png");

    // tag
    // efficiency using different tag-pt cut: data and MC
    CPlot plotEffEtaPt_tagPosDat((outDir+"/plots_"+eType+"/effetapt_tagPosDat").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_tagPosDat.SetLogy(doLogy);
    plotEffEtaPt_tagPosDat.AddHist2D(tagPosDat, "COLZ,text");
    plotEffEtaPt_tagPosDat.Draw(c, 1, "png");
    CPlot plotEffEtaPt_tagPosMC((outDir+"/plots_"+eType+"/effetapt_tagPosMC").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_tagPosMC.SetLogy(doLogy);
    plotEffEtaPt_tagPosMC.AddHist2D(tagPosMC, "COLZ,text");
    plotEffEtaPt_tagPosMC.Draw(c, 1, "png");
    // Data ratio between two tag pt cuts
    tagPosDat->Divide(amcPosDat);
    CPlot plotEffEtaPt_tagPosDatRatio((outDir+"/plots_"+eType+"/effetapt_tagPosDatRatio").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_tagPosDatRatio.SetLogy(doLogy);
    plotEffEtaPt_tagPosDatRatio.AddHist2D(tagPosDat, "COLZ,text");
    plotEffEtaPt_tagPosDatRatio.Draw(c, 1, "png");
    // MC ratio between two tag pt cuts
    tagPosMC->Divide(amcPosMC);
    CPlot plotEffEtaPt_tagPosMCRatio((outDir+"/plots_"+eType+"/effetapt_tagPosMCRatio").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_tagPosMCRatio.SetLogy(doLogy);
    plotEffEtaPt_tagPosMCRatio.AddHist2D(tagPosMC, "COLZ,text");
    plotEffEtaPt_tagPosMCRatio.Draw(c, 1, "png");
    // double ratio betwee data fits and MC fits
    tagPosDat->Divide(tagPosMC);
    CPlot plotEffEtaPt_tagPosDatMCRatio((outDir+"/plots_"+eType+"/effetapt_tagPosDatMCRatio").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_tagPosDatMCRatio.SetLogy(doLogy);
    plotEffEtaPt_tagPosDatMCRatio.AddHist2D(tagPosDat, "COLZ,text");
    plotEffEtaPt_tagPosDatMCRatio.Draw(c, 1, "png");

    // data ratio between of the toy MC using alternative model and default model
    hEffTagPos->Divide(amcPosDat);
    CPlot plotEffEtaPt_tagPosToyDatRatio((outDir+"/plots_"+eType+"/effetapt_tagPosToyDatRatio").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_tagPosToyDatRatio.SetLogy(doLogy);
    plotEffEtaPt_tagPosToyDatRatio.AddHist2D(hEffTagPos, "COLZ,text");
    plotEffEtaPt_tagPosToyDatRatio.Draw(c, 1, "png");

    // double ratios between toyMC on data and MC
    CPlot plotEffEtaPt_tagPosSys((outDir+"/plots_"+eType+"/effetapt_tagPosSys").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_tagPosSys.SetLogy(doLogy);
    plotEffEtaPt_tagPosSys.AddHist2D(hTagPos, "COLZ,text");
    plotEffEtaPt_tagPosSys.Draw(c, 1, "png");

    // original ratios
    amcPosDat->Divide(amcPosMC);
    CPlot plotEffEtaPt_amcPosRatio((outDir+"/plots_"+eType+"/effetapt_amcCorr").Data(), "", "probe #eta", "probe p_{T} [GeV/c]");
    plotEffEtaPt_amcPosRatio.SetLogy(doLogy);
    plotEffEtaPt_amcPosRatio.AddHist2D(amcPosDat, "COLZ,text");
    plotEffEtaPt_amcPosRatio.Draw(c, 1, "png");

    makeHTML(eType);

    char outfilename[500];
    sprintf(outfilename, "%s/SysUnc_%s.root", outDir.Data(), eType.Data());
    TFile* f = new TFile(outfilename, "RECREATE");
    hFSRNeg->Write();
    hFSRPos->Write();
    hMCNeg->Write();
    hMCPos->Write();
    hBkgNeg->Write();
    hBkgPos->Write();
    hTagNeg->Write();
    hTagPos->Write();
    f->Close();
}
