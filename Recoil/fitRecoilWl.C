//================================================================================================
//
// Perform fits to recoil against W->lnu events
//
//  * Outputs a ROOT file of functions parametrizing the distribution of recoil components
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TLorentzVector.h" // 4-vector class
#include <TF1.h>            // 1D function
#include <TFile.h>          // file handle class
#include <TFitResult.h>     // class to handle fit results
#include <TGraphErrors.h>   // graph class
#include <TTree.h>          // class to access ntuples
#include <fstream>          // standard I/O
#include <iostream>         // standard I/O
#include <sstream>

#include "MitEwk13TeV/Utils/CPlot.hh"          // helper class for plots
#include "MitEwk13TeV/Utils/MitStyleRemix.hh"  // style settings for drawing
#include "MitEwk13TeV/Utils/METXYCorrector.hh" // MET XY correction
#include "MitEwk13TeV/Utils/MyTools.hh"        // various helper functions
#include "MitEwk13TeV/RochesterCorr/RoccoR.cc"
#include "MitEwk13TeV/Utils/AppEffSF.cc"

#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooRealIntegral.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TRandom.h"
#endif

using namespace RooFit;
using namespace std;

bool doElectron = false;

//=== FUNCTION DECLARATIONS ======================================================================================
//--------------------------------------------------------------------------------------------------
// perform fit of recoil component
void performFit(const vector<TH1D *> hv, const vector<TH1D *> hbkgv, const Double_t *ptbins, const Int_t nbins,
                const Int_t model, const Bool_t sigOnly,
                const vector<RooDataSet> lDataSet, const vector<RooRealVar> lVar,
                TCanvas *c, const char *plabel, const char *xlabel,
                Double_t *mean1Arr, Double_t *mean1ErrArr,
                Double_t *mean2Arr, Double_t *mean2ErrArr,
                Double_t *mean3Arr, Double_t *mean3ErrArr,
                Double_t *sigma0Arr, Double_t *sigma0ErrArr,
                Double_t *sigma1Arr, Double_t *sigma1ErrArr,
                Double_t *sigma2Arr, Double_t *sigma2ErrArr,
                Double_t *sigma3Arr, Double_t *sigma3ErrArr,
                Double_t *frac2Arr, Double_t *frac2ErrArr,
                Double_t *frac3Arr, Double_t *frac3ErrArr,
                RooWorkspace *workspace,
                const char *outputDir,
                int etaBinCategory, bool do_keys,
                bool do_5TeV);

//=== MAIN MACRO =================================================================================================

void fitRecoilWl(TString indir,                    // input ntuple
                 Int_t pfu1model,                  // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
                 Int_t pfu2model,                  // u2 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
                 Bool_t sigOnly,                   // signal event only?
                 Int_t charge,                     // charge requirement
                 Bool_t useData,                   // use Data? (0 = signal MC, 1 = data)
                 std::string metVar = "met",       // variable name to pull
                 std::string metPhiVar = "metPhi", // variable name for MET phi
                 std::string metName = "pf",       // for recordkeeping?
                 TString outputDir = "./",         // output directory
                 Double_t lumi = 1,                // lumi value for the data fits
                 int etaBinCategory = 0,           // 0 is inclusive, 1 is fabs(eta)<=0.5,  2 is fabs(eta)=[0.5,1], 3 is fabs(eta)>=1
                 bool do_keys = 0,                 // 0 for regular function, 1 for RooKeysPDF
                 bool do_5TeV = 0,                 // 0 for 13 TeV, 1 for 5 TeV
                 bool doElectron = 0               // 0 for muon, 1 for electron
)
{

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================
    // preserving the fine binning at low pT but the higher-pT bins (>75 GeV have been adjusted to be slightly wider)
    Double_t ptbins[] = {0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 65, 70, 75, 80, 90, 100, 125, 150, 1000};
    // Double_t ptbins[] = {0, 25, 50, 75, 100, 125, 150, 1000};
    Int_t nbins = sizeof(ptbins) / sizeof(Double_t) - 1;

    vector<TString> fnamev;
    vector<Bool_t> isBkgv;

    if (useData == 0)
    {
        if (doElectron and !do_5TeV)
        {

            fnamev.push_back(indir + "/we0.root");
            isBkgv.push_back(kFALSE);
            fnamev.push_back(indir + "/we1.root");
            isBkgv.push_back(kFALSE);
            fnamev.push_back(indir + "/we2.root");
            isBkgv.push_back(kFALSE);
        }
        else if (doElectron and do_5TeV)
        {
            fnamev.push_back(indir + "/we.root");
            isBkgv.push_back(kFALSE);
        }
        else if (!doElectron and !do_5TeV)
        {
            fnamev.push_back(indir + "/wm0.root");
            isBkgv.push_back(kFALSE);
            fnamev.push_back(indir + "/wm1.root");
            isBkgv.push_back(kFALSE);
            fnamev.push_back(indir + "/wm2.root");
            isBkgv.push_back(kFALSE);
        }
        else
        {
            // 5TeV muon
            fnamev.push_back(indir + "/wm.root");
            isBkgv.push_back(kFALSE);
        }
    }
    else
    {
        fnamev.push_back(TString(indir) + TString("/data.root"));
        isBkgv.push_back(kFALSE);
    }

    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;
    Double_t mu_MASS = 0.1057;
    if (doElectron)
    {
        // call electrons "muons"
        mu_MASS = 0.000511;
    }
    const Double_t ETA_ECAL_GAP_LOW = 1.4442;
    const Double_t ETA_ECAL_GAP_HIGH = 1.566;

    const bool doMETXYCorrection = true;

    // Setting up rochester corrections
    const TString envStr = (TString)gSystem->Getenv("CMSSW_BASE") + "/src/";
    RoccoR rc((envStr + "/MitEwk13TeV/RochesterCorr/RoccoR2017.txt").Data());

    TString sqrts = "13TeV";
    if (do_5TeV)
        sqrts = "5TeV";
    TString channel = "mumu";
    if (doElectron)
        channel = "ee";
    METXYCorrector metXYCorr("XYCorrector", (envStr + "/MitEwk13TeV/Recoil/data/met_xy_" + sqrts + "_" + channel + ".root").Data());

    // load efficiency corrections
    channel = "Zmm";
    if (doElectron)
        channel = "Zee";
    AppEffSF effs("");
    TString baseDir = envStr + "Corrections/Efficiency/" + sqrts + "/results/" + channel + "/";
    effs.SetBaseDir(baseDir);
    if (doElectron)
    {
        effs.loadHLT("EleHLTEff_aMCxPythia", "Positive", "Negative");
        effs.loadSel("EleGSFSelEff_aMCxPythia", "Combined", "Combined");
    }
    else
    {
        effs.loadHLT("MuHLTEff_aMCxPythia", "Positive", "Negative");
        effs.loadSel("MuSITEff_aMCxPythia", "Combined", "Combined");
        effs.loadSta("MuStaEff_aMCxPythia", "Combined", "Combined");
    }

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    char hname[100];
    vector<TH1D *> hPFu1v, hPFu1Bkgv;
    vector<TH1D *> hPFu2v, hPFu2Bkgv;

    vector<RooDataSet> lDataSetU1;
    vector<RooDataSet> lDataSetU2;

    vector<RooRealVar> vu1Var;
    vector<RooRealVar> vu2Var;

    RooWorkspace pdfsU1("pdfsU1");
    RooWorkspace pdfsU2("pdfsU2");

    for (Int_t ibin = 0; ibin < nbins; ibin++)
    {

        int range = 60;
        if (ptbins[ibin] > 50)
            range = 80;
        if (ptbins[ibin] > 80)
            range = 100;
        if (ptbins[ibin] > 120)
            range = 150;
        sprintf(hname, "hPFu1_%i", ibin);
        // hPFu1v.push_back(new TH1D(hname, "", 100, -range - ptbins[ibin], range - ptbins[ibin]));
        hPFu1v.push_back(new TH1D(hname, "", 100, -range, range));
        hPFu1v[ibin]->Sumw2();
        sprintf(hname, "hPFu1Bkg_%i", ibin);
        // hPFu1Bkgv.push_back(new TH1D(hname, "", 100, -range - ptbins[ibin], range - ptbins[ibin]));
        hPFu1Bkgv.push_back(new TH1D(hname, "", 100, -range, range));
        hPFu1Bkgv[ibin]->Sumw2();

        sprintf(hname, "hPFu2_%i", ibin);
        hPFu2v.push_back(new TH1D(hname, "", 100, -range, range));
        hPFu2v[ibin]->Sumw2();
        sprintf(hname, "hPFu2Bkg_%i", ibin);
        hPFu2Bkgv.push_back(new TH1D(hname, "", 100, -range, range));
        hPFu2Bkgv[ibin]->Sumw2();

        std::stringstream name;
        name << "u_" << ibin;
        // RooRealVar u1Var(name.str().c_str(), name.str().c_str(), 0, -range - ptbins[ibin], range - ptbins[ibin]);
        RooRealVar u1Var(name.str().c_str(), name.str().c_str(), 0, -range, range);
        RooRealVar u2Var(name.str().c_str(), name.str().c_str(), 0, -range, range);

        name.str("");
        name << "weight_" << ibin;
        // weightVar
        RooRealVar weightVar(name.str().c_str(), name.str().c_str(), -10000, 10000);

        vu1Var.push_back(u1Var);
        vu2Var.push_back(u2Var);

        sprintf(hname, "hDataSetU1_%i", ibin);
        RooDataSet dataSetU1(hname, hname, RooArgSet(u1Var, weightVar), RooFit::WeightVar(weightVar));
        std::cout << " is Weighted " << dataSetU1.isWeighted() << std::endl;
        lDataSetU1.push_back(dataSetU1);
        sprintf(hname, "hDataSetU2_%i", ibin);
        RooDataSet dataSetU2(hname, hname, RooArgSet(u2Var, weightVar), RooFit::WeightVar(weightVar));
        lDataSetU2.push_back(dataSetU2);
    }

    TFile *infile = 0;
    TTree *intree = 0;

    //
    // Declare output ntuple variables
    //
    UInt_t runNum, lumiSec, evtNum;
    // UInt_t  npv, npu;
    Float_t scale1fb, puWeight; //, scale1fbUp, scale1fbDown;
    Float_t met, metPhi;        //, mt, u1, u2;
    Int_t q;
    Float_t prefireWeight;
    UInt_t nTkLayers = 0;  // for roch corr
    Float_t genMuonPt = 0; // gen muon pt for roch corr
    TLorentzVector *lep = 0, *lep_raw = 0, *genV = 0, *genLep = 0;

    //   Float_t puWeight;
    //   Float_t scale1fb;

    for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++)
    {
        cout << "Processing " << fnamev[ifile] << "..." << endl;
        bool isData = (fnamev[ifile].Contains("data.root"));
        cout << "IsData = " << isData << endl;

        infile = new TFile(fnamev[ifile]);
        intree = (TTree *)infile->Get("Events");

        intree->SetBranchAddress("scale1fb", &scale1fb);           // event weight per 1/fb (MC)
        intree->SetBranchAddress(metVar.c_str(), &met);            // Uncorrected PF MET
        intree->SetBranchAddress(metPhiVar.c_str(), &metPhi);      // phi(MET)
        intree->SetBranchAddress("q", &q);                         // lepton charge
        intree->SetBranchAddress("lep", &lep);                     // lepton 4-vector
        intree->SetBranchAddress("prefireWeight", &prefireWeight); // prefire weight for 2017 conditions (MC)
        if (!doElectron)
        {
            intree->SetBranchAddress("nTkLayers", &nTkLayers); // lepton 4-vector
            intree->SetBranchAddress("genMuonPt", &genMuonPt); // GEN muon pt (signal MC
        }
        intree->SetBranchAddress("genV", &genV);     // GEN W boson 4-vector (signal MC)
        intree->SetBranchAddress("genLep", &genLep); // GEN lepton 4-vector (signal MC)
        if (doElectron)
            intree->SetBranchAddress("lep_raw", &lep_raw); // probe lepton 4-vector

        TH1D *hGenWeights;
        double totalNorm = 1.0;
        if (!isData)
        {
            hGenWeights = (TH1D *)infile->Get("hGenWeights");
            totalNorm = hGenWeights->Integral();
        }

        int iterator = 1;
        if (do_keys and sigOnly and !do_5TeV)
        {
            // to speed up the RooKeysPdf, we only use 1/iterator of the events
            // otherwise it takes too long
            // Change the totalNorm accordingly
            if (ifile == 0)
                iterator = 20;
            else if (ifile == 1)
                iterator = 5;
            else
                iterator = 2;
        }
        // if (do_5TeV and sigOnly)
        //     iterator = 2;
        totalNorm = totalNorm / iterator;

        //
        // Loop over events
        //
        for (Int_t ientry = 0; ientry < intree->GetEntries(); ientry += iterator)
        {
            intree->GetEntry(ientry);

            if (ientry % 10000 == 0)
                cout << "Processing event " << ientry << endl;

            // apply rochester correction
            TLorentzVector mu;
            mu.SetPtEtaPhiM(lep->Pt(), lep->Eta(), lep->Phi(), mu_MASS);
            double SF1 = 1;
            if (doElectron)
            {
                if (isData)
                {
                    SF1 = rc.kScaleDT(q, mu.Pt(), mu.Eta(), mu.Phi());
                }
                else
                {
                    if (genMuonPt > 0)
                        SF1 = rc.kSpreadMC(q, mu.Pt(), mu.Eta(), mu.Phi(), genMuonPt);
                    else
                        SF1 = rc.kSmearMC(q, mu.Pt(), mu.Eta(), mu.Phi(), nTkLayers, gRandom->Uniform(1));
                }
                mu *= SF1;
            }

            if (charge == 1 && q < 0)
                continue;
            if (charge == -1 && q > 0)
                continue;
            if (mu.Pt() < PT_CUT)
                continue;
            if (fabs(lep->Eta()) > ETA_CUT)
                continue;
            if (doElectron)
            {
                if (fabs(lep->Eta()) >= ETA_ECAL_GAP_LOW && fabs(lep->Eta()) <= ETA_ECAL_GAP_HIGH)
                    continue;
            }

            float genVy = genV->Rapidity();
            float genVPt = genV->Pt();
            float genVPhi = genV->Phi();

            // 0 is inclusive, 1 is fabs(eta)<=0.5,  2 is fabs(eta)=[0.5,1], 3 is fabs(eta)>=1
            if (etaBinCategory == 1 && fabs(genVy) > 0.5)
                continue;
            if (etaBinCategory == 2 && (fabs(genVy) <= 0.5 || fabs(genVy) >= 1))
                continue;
            if (etaBinCategory == 3 && fabs(genVy) < 1)
                continue;

            Int_t ipt = -1;
            for (Int_t ibin = 0; ibin < nbins; ibin++)
            {
                if (genVPt > ptbins[ibin] && genVPt <= ptbins[ibin + 1])
                    ipt = ibin;
            }
            if (ipt < 0)
                continue;

            // run the MET XY Correction first
            if (doMETXYCorrection)
            {
                // std::cout << "before correction met = " << met << " metPhi = " << metPhi << std::endl;
                metXYCorr.CorrectMETXY(met, metPhi, isData);
                // std::cout << "after correction met = " << met << " metPhi = " << metPhi << std::endl;
            }

            TVector2 vLepRaw1;
            if (doElectron)
            {
                vLepRaw1.Set((lep_raw->Pt()) * cos(lep_raw->Phi()), (lep_raw->Pt()) * sin(lep_raw->Phi()));
            }
            else
            {
                vLepRaw1.Set((lep->Pt()) * cos(lep->Phi()), (lep->Pt()) * sin(lep->Phi()));
            }
            TVector2 vLepCor1((mu.Pt()) * cos(lep->Phi()), (mu.Pt()) * sin(lep->Phi()));
            TVector2 vMetCorr((met)*cos(metPhi), (met)*sin(metPhi));
            Double_t corrMetWithLepton = (vMetCorr + vLepRaw1 - vLepCor1).Mod();
            Double_t corrMetWithLeptonPhi = (vMetCorr + vLepRaw1 - vLepCor1).Phi();
            double pUX = corrMetWithLepton * cos(corrMetWithLeptonPhi) + mu.Pt() * cos(lep->Phi());
            double pUY = corrMetWithLepton * sin(corrMetWithLeptonPhi) + mu.Pt() * sin(lep->Phi());
            double pU = sqrt(pUX * pUX + pUY * pUY);
            // projected on the corrected W
            double pCos = -(pUX * cos(genVPhi) + pUY * sin(genVPhi)) / pU;
            double pSin = (pUX * sin(genVPhi) - pUY * cos(genVPhi)) / pU;

            double pU1 = pU * pCos; // U1 in data
            double pU2 = pU * pSin; // U2 in data
            pU1 = pU1 + genVPt;

            vu1Var[ipt].setVal(pU1);
            vu2Var[ipt].setVal(pU2);

            double corr = 1.0;
            if (!isData)
                corr = effs.fullCorrections(&mu, 1);

            double weight = 1.0;
            if (!isData)
                weight = scale1fb * lumi * corr * prefireWeight / totalNorm;

            if (isBkgv[ifile] && !sigOnly)
            {
                // bkg contribute negatively
                weight = weight * (-1);
            }
            else if (isBkgv[ifile] && sigOnly)
            {
                weight = 0;
            }
            lDataSetU1[ipt].add(RooArgSet(vu1Var[ipt]), weight); // need to add the weights
            lDataSetU2[ipt].add(RooArgSet(vu2Var[ipt]), weight);

            if (ientry < 200)
                std::cout << "weight = " << weight << " scale1fb " << scale1fb << " lumi " << lumi << " totalNorm " << totalNorm << std::endl;

            if (isBkgv[ifile])
            {
                hPFu1Bkgv[ipt]->Fill(pU1, scale1fb * lumi / totalNorm);
                hPFu2Bkgv[ipt]->Fill(pU2, scale1fb * lumi / totalNorm);
            }
            else
            {
                hPFu1v[ipt]->Fill(pU1, scale1fb * lumi / totalNorm);
                hPFu2v[ipt]->Fill(pU2, scale1fb * lumi / totalNorm);
            }
        }

        delete infile;
        infile = 0, intree = 0;
    }

    Double_t pfu1Mean[nbins], pfu1MeanErr[nbins];
    Double_t pfu1Mean2[nbins], pfu1Mean2Err[nbins];
    Double_t pfu1Mean3[nbins], pfu1Mean3Err[nbins];
    Double_t pfu1Sigma0[nbins], pfu1Sigma0Err[nbins];
    Double_t pfu1Sigma1[nbins], pfu1Sigma1Err[nbins];
    Double_t pfu1Sigma2[nbins], pfu1Sigma2Err[nbins];
    Double_t pfu1Sigma3[nbins], pfu1Sigma3Err[nbins];
    Double_t pfu1Frac2[nbins], pfu1Frac2Err[nbins];
    Double_t pfu1Frac3[nbins], pfu1Frac3Err[nbins];

    Double_t pfu2Mean[nbins], pfu2MeanErr[nbins];
    Double_t pfu2Mean2[nbins], pfu2Mean2Err[nbins];
    Double_t pfu2Mean3[nbins], pfu2Mean3Err[nbins];
    Double_t pfu2Sigma0[nbins], pfu2Sigma0Err[nbins];
    Double_t pfu2Sigma1[nbins], pfu2Sigma1Err[nbins];
    Double_t pfu2Sigma2[nbins], pfu2Sigma2Err[nbins];
    Double_t pfu2Sigma3[nbins], pfu2Sigma3Err[nbins];
    Double_t pfu2Frac2[nbins], pfu2Frac2Err[nbins];
    Double_t pfu2Frac3[nbins], pfu2Frac3Err[nbins];

    TCanvas *c = MakeCanvas("c", "c", 800, 800);

    // Fitting PF-MET u1
    performFit(hPFu1v, hPFu1Bkgv, ptbins, nbins, pfu1model, sigOnly,
               lDataSetU1, vu1Var,
               c, "pfu1", "u_{#parallel  } + p^{W}_{T} [GeV]",
               pfu1Mean, pfu1MeanErr,
               pfu1Mean2, pfu1Mean2Err,
               pfu1Mean3, pfu1Mean3Err,
               pfu1Sigma0, pfu1Sigma0Err,
               pfu1Sigma1, pfu1Sigma1Err,
               pfu1Sigma2, pfu1Sigma2Err,
               pfu1Sigma3, pfu1Sigma3Err,
               pfu1Frac2, pfu1Frac2Err,
               pfu1Frac3, pfu1Frac3Err,
               &pdfsU1,
               outputDir,
               etaBinCategory, do_keys,
               do_5TeV);
    char outpdfname[50];
    sprintf(outpdfname, "%s/%s.root", outputDir.Data(), "pdfsU1");
    pdfsU1.writeToFile(outpdfname);

    // Fitting PF-MET u2
    performFit(hPFu2v, hPFu2Bkgv, ptbins, nbins, pfu2model, sigOnly,
               lDataSetU2, vu2Var,
               c, "pfu2", "u_{#perp  } [GeV]",
               pfu2Mean, pfu2MeanErr,
               pfu2Mean2, pfu2Mean2Err,
               pfu2Mean3, pfu2Mean3Err,
               pfu2Sigma0, pfu2Sigma0Err,
               pfu2Sigma1, pfu2Sigma1Err,
               pfu2Sigma2, pfu2Sigma2Err,
               pfu2Sigma3, pfu2Sigma3Err,
               pfu2Frac2, pfu2Frac2Err,
               pfu2Frac3, pfu2Frac3Err,
               &pdfsU2,
               outputDir,
               etaBinCategory, do_keys,
               do_5TeV);
    sprintf(outpdfname, "%s/%s.root", outputDir.Data(), "pdfsU2");
    pdfsU2.writeToFile(outpdfname);

    // Save histograms
    char outhistname[50];
    sprintf(outhistname, "%s/%s.root", outputDir.Data(), "histos");
    TFile *ofile = new TFile(outhistname, "recreate");
    for (Int_t ibin = 0; ibin < nbins; ibin++)
    {
        hPFu1v[ibin]->Write();
        hPFu2v[ibin]->Write();
        hPFu1Bkgv[ibin]->Write();
        hPFu2Bkgv[ibin]->Write();
    }
    ofile->Close();

    delete infile;
    infile = 0, intree = 0;

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;
}

//--------------------------------------------------------------------------------------------------
void performFit(const vector<TH1D *> hv, const vector<TH1D *> hbkgv, const Double_t *ptbins, const Int_t nbins,
                const Int_t model, const Bool_t sigOnly,
                vector<RooDataSet> lDataSet, vector<RooRealVar> lVar,
                TCanvas *c, const char *plabel, const char *xlabel,
                Double_t *mean1Arr, Double_t *mean1ErrArr,
                Double_t *mean2Arr, Double_t *mean2ErrArr,
                Double_t *mean3Arr, Double_t *mean3ErrArr,
                Double_t *sigma0Arr, Double_t *sigma0ErrArr,
                Double_t *sigma1Arr, Double_t *sigma1ErrArr,
                Double_t *sigma2Arr, Double_t *sigma2ErrArr,
                Double_t *sigma3Arr, Double_t *sigma3ErrArr,
                Double_t *frac2Arr, Double_t *frac2ErrArr,
                Double_t *frac3Arr, Double_t *frac3ErrArr,
                RooWorkspace *wksp,
                const char *outputDir,
                int etaBinCategory, bool do_keys,
                bool do_5TeV)
{
    char pname[50];
    char lumi[200];
    char ylabel[50];
    char binlabel[50];
    char binYlabel[50];
    char nsigtext[50];
    char nbkgtext[50];
    char chi2text[50];
    char hmeantext[50];
    char hrmstext[50];

    char mean1text[50];
    char mean2text[50];
    char mean3text[50];
    char sig0text[50];
    char sig1text[50];
    char sig2text[50];
    char sig3text[50];
    char frac2text[50];
    char frac3text[50];

    for (Int_t ibin = 0; ibin < nbins; ibin++)
    {

        std::stringstream name;
        // unfortunately have to give each variable individual names for each bin
        name << "u_" << ibin;
        RooRealVar u(name.str().c_str(), name.str().c_str(), hv[ibin]->GetXaxis()->GetXmin(), hv[ibin]->GetXaxis()->GetXmax());
        name.str("");
        u.setBins(100);
        RooDataHist dataHist("dataHist", "dataHist", RooArgSet(u), hv[ibin]);

        //
        // Set up background histogram templates
        //
        RooDataHist bkgHist("bkgHist", "bkgHist", RooArgSet(u), hbkgv[ibin]);
        RooHistPdf bkg("bkg", "bkg", u, bkgHist, 0);

        //
        // Set up fit parameters
        //
        name.str("");
        name << "mean1_" << ibin;
        RooRealVar mean1(name.str().c_str(), name.str().c_str(),
                         hv[ibin]->GetMean(),
                         hv[ibin]->GetXaxis()->GetXmin() + 50,
                         hv[ibin]->GetXaxis()->GetXmax() - 50);
        name.str("");
        name << "mean2_" << ibin;
        RooRealVar mean2(name.str().c_str(), name.str().c_str(),
                         hv[ibin]->GetMean(),
                         hv[ibin]->GetXaxis()->GetXmin() + 50,
                         hv[ibin]->GetXaxis()->GetXmax() - 50);
        name.str("");
        name << "mean3_" << ibin;
        RooRealVar mean3(name.str().c_str(), name.str().c_str(),
                         hv[ibin]->GetMean() * 0.85,
                         hv[ibin]->GetXaxis()->GetXmin() + 50,
                         hv[ibin]->GetXaxis()->GetXmax() - 50);
        name.str("");
        name << "sigma1_" << ibin;
        RooRealVar sigma1(name.str().c_str(), name.str().c_str(), 0.3 * (hv[ibin]->GetRMS()), 0., 2.3 * (hv[ibin]->GetRMS()));
        name.str("");
        name << "sigma2_" << ibin;
        RooRealVar sigma2(name.str().c_str(), name.str().c_str(), 1.0 * (hv[ibin]->GetRMS()), 0., 4.5 * (hv[ibin]->GetRMS()));
        name.str("");
        name << "sigma3_" << ibin;
        RooRealVar sigma3(name.str().c_str(), name.str().c_str(), 2.0 * (hv[ibin]->GetRMS()), 0., 9 * (hv[ibin]->GetRMS()));
        name.str("");
        name << "frac2_" << ibin;
        RooRealVar frac2(name.str().c_str(), name.str().c_str(), 0.5, 0.15, 0.85);
        name.str("");
        name << "frac3_" << ibin;
        RooRealVar frac3(name.str().c_str(), name.str().c_str(), 0.05, 0, 0.15);

        if (model == 2)
        {
            frac2.setVal(0.5);
            frac2.setRange(0., 1.);
        }

        if (string(plabel) == string("pfu2"))
        {
            mean1.setVal(0);
            mean1.setConstant(kTRUE);
            mean2.setVal(0);
            mean2.setConstant(kTRUE);
            mean3.setVal(0);
            mean3.setConstant(kTRUE);
        }

        name.str("");
        name << "gauss1_" << ibin;
        RooGaussian gauss1(name.str().c_str(), name.str().c_str(), u, mean1, sigma1);
        name.str("");
        name << "gauss2_" << ibin;
        // RooGaussian gauss2(name.str().c_str(), name.str().c_str(), u, mean2, sigma2);
        RooGaussian gauss2(name.str().c_str(), name.str().c_str(), u, mean1, sigma2);
        name.str("");
        name << "gauss3_" << ibin;
        RooGaussian gauss3(name.str().c_str(), name.str().c_str(), u, mean3, sigma3);
        name.str("");

        RooGaussian constGauss1("constGauss1", "constGauss1", mean1, RooConst(hv[ibin]->GetMean()), RooConst(0.15 * hv[ibin]->GetRMS()));
        RooGaussian constGauss2("constGauss2", "constGauss2", mean2, RooConst(hv[ibin]->GetMean()), RooConst(0.15 * hv[ibin]->GetRMS()));
        RooGaussian constGauss3("constGauss3", "constGauss3", mean3, RooConst(hv[ibin]->GetMean()), RooConst(0.15 * hv[ibin]->GetRMS()));

        //
        // Define formula for overall width (sigma0)
        //
        char formula[100];
        RooArgList params;
        if (model == 1)
        {
            sprintf(formula, "sigma1");
        }
        else if (model == 2)
        {
            sprintf(formula, "(1.-frac2)*sigma1 + frac2*sigma2");
            params.add(frac2);
            params.add(sigma1);
            params.add(sigma2);
        }
        else if (model == 3)
        {
            sprintf(formula, "(1.-frac2-frac3)*sigma1 + frac2*sigma2 + frac3*sigma3");
            params.add(frac2);
            params.add(frac3);
            params.add(sigma1);
            params.add(sigma2);
            params.add(sigma3);
        }
        RooFormulaVar sigma0("sigma0", formula, params);

        //
        // Construct fit model
        //
        RooArgList shapes;
        if (model >= 3)
            shapes.add(gauss3);
        if (model >= 2)
            shapes.add(gauss2);
        shapes.add(gauss1);

        RooArgList fracs;
        if (model >= 3)
            fracs.add(frac3);
        if (model >= 2)
            fracs.add(frac2);

        // sig pdfsU1
        name.str("");
        name << "sig_" << ibin;
        // recursive: sig = frac2 * gauss2 + (1-frac2) * gauss1
        RooAddPdf sig(name.str().c_str(), name.str().c_str(), shapes, fracs, true);
        name.str("");

        RooArgList parts;
        parts.add(sig);
        if (!sigOnly)
            parts.add(bkg);

        RooArgList yields;
        name.str("");
        name << "nsig_" << ibin;
        RooRealVar nsig(name.str().c_str(), name.str().c_str(), 0.98 * (hv[ibin]->Integral()), 0., 1.1 * hv[ibin]->Integral()); // just to be sure that doesn't it the boundary
        yields.add(nsig);
        name.str("");
        name << "nbkg_" << ibin;
        RooRealVar nbkg(name.str().c_str(), name.str().c_str(), 0.01 * (hv[ibin]->Integral()), 0., 0.50 * (hv[ibin]->Integral()));
        if (!sigOnly)
            yields.add(nbkg);
        else
            nbkg.setVal(0);

        //     std::stringstream name;
        name.str("");
        name << "modelpdf_" << ibin << std::endl;
        RooAddPdf modelpdf(name.str().c_str(), name.str().c_str(), parts, yields);
        name.str("");

        std::cout << "name = " << name.str().c_str() << std::endl;

        std::cout << "on bin # " << ibin << std::endl;
        std::cout << "signal evts = " << nsig.getVal() << std::endl;
        std::cout << "total events = " << hv[ibin]->Integral() << std::endl;
        std::cout << "sigma1 = " << sigma1.getVal() << std::endl;
        std::cout << "sigma2 = " << sigma2.getVal() << std::endl;

        if (ibin > 0)
            std::cout << "sigma max = " << 1.5 * hv[ibin - 1]->GetRMS() << " " << 1.8 * hv[ibin - 1]->GetRMS() << std::endl;

        //
        // Perform fit
        //
        RooFitResult *fitResult = 0;
        fitResult = modelpdf.fitTo(dataHist,
                                   NumCPU(4),
                                   SumW2Error(kTRUE),
                                   Minimizer("Minuit2", "minimize"),
                                   ExternalConstraints(constGauss1), ExternalConstraints(constGauss2),
                                   // ExternalConstraints(constGauss3),
                                   RooFit::Minos(),
                                   RooFit::Strategy(2),
                                   RooFit::Save());

        int nTries = 0;
        do
        {
            // if(ibin==22||ibin==27)break;
            fitResult = modelpdf.fitTo(dataHist,
                                       NumCPU(4),
                                       SumW2Error(kTRUE),
                                       Minimizer("Minuit2", "scan"),
                                       ExternalConstraints(constGauss1), ExternalConstraints(constGauss2),
                                       // ExternalConstraints(constGauss3),
                                       RooFit::Minos(),
                                       RooFit::Strategy(2),
                                       RooFit::Save());
            fitResult = modelpdf.fitTo(dataHist,
                                       NumCPU(4),
                                       SumW2Error(kTRUE),
                                       Minimizer("Minuit2", "migrad"),
                                       ExternalConstraints(constGauss1), ExternalConstraints(constGauss2),
                                       // ExternalConstraints(constGauss3),
                                       RooFit::Hesse(),
                                       RooFit::Strategy(2),
                                       RooFit::Save());
            fitResult = modelpdf.fitTo(dataHist,
                                       NumCPU(4),
                                       SumW2Error(kTRUE),
                                       Minimizer("Minuit2", "improve"),
                                       ExternalConstraints(constGauss1), ExternalConstraints(constGauss2),
                                       // ExternalConstraints(constGauss3),
                                       RooFit::Minos(),
                                       RooFit::Strategy(2),
                                       RooFit::Save());
            fitResult = modelpdf.fitTo(dataHist,
                                       NumCPU(4),
                                       SumW2Error(kTRUE),
                                       Minimizer("Minuit2", "minimize"),
                                       ExternalConstraints(constGauss1), ExternalConstraints(constGauss2),
                                       // ExternalConstraints(constGauss3),
                                       RooFit::Minos(),
                                       RooFit::Strategy(2),
                                       RooFit::Save());
            nTries++;
        } while ((fitResult->status() > 0 || fitResult->covQual() < 3) && nTries < 2);

        c->SetFillColor(kWhite);
        if (fitResult->status() > 0)
            c->SetFillColor(kYellow);

        wksp->import(u);
        wksp->import(modelpdf);

        mean1Arr[ibin] = mean1.getVal();
        mean1ErrArr[ibin] = mean1.getError();
        sigma1Arr[ibin] = sigma1.getVal();
        sigma1ErrArr[ibin] = sigma1.getError();
        if (model >= 2)
        {
            mean2Arr[ibin] = mean2.getVal();
            mean2ErrArr[ibin] = mean2.getError();
            sigma0Arr[ibin] = sigma0.getVal();
            sigma0ErrArr[ibin] = sigma0.getPropagatedError(*fitResult);
            sigma2Arr[ibin] = sigma2.getVal();
            sigma2ErrArr[ibin] = sigma2.getError();
            frac2Arr[ibin] = frac2.getVal();
            frac2ErrArr[ibin] = frac2.getError();
        }
        if (model >= 3)
        {
            mean3Arr[ibin] = mean3.getVal();
            mean3ErrArr[ibin] = mean3.getError();
            sigma3Arr[ibin] = sigma3.getVal();
            sigma3ErrArr[ibin] = sigma3.getError();
            frac3Arr[ibin] = frac3.getVal();
            frac3ErrArr[ibin] = frac3.getError();
        }

        std::cout << "Plot Fit results " << std::endl;
        //
        // Plot fit results
        //
        RooPlot *frame = u.frame(Bins(100));
        dataHist.plotOn(frame, MarkerStyle(kFullCircle), MarkerSize(0.8), DrawOption("ZP"));
        modelpdf.plotOn(frame);
        frame->GetYaxis()->SetTitleOffset(1.2);

        if (!sigOnly)
            modelpdf.plotOn(frame, Components("bkg"), LineStyle(kDotted), LineColor(kMagenta + 2));
        name.str("");
        name << "gauss1_" << ibin;
        if (model >= 2)
            sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kRed));
        name.str("");
        name << "gauss2_" << ibin;
        if (model >= 2)
            sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kMagenta));
        name.str("");
        name << "gauss3_" << ibin;
        if (model >= 3)
            sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kGreen + 2));

        // draw the curve
        sig.plotOn(frame, FillColor(7), VisualizeError(*fitResult, 1), RooFit::Components(sig)); // 1 sigma band
        sig.plotOn(frame, RooFit::LineColor(kBlue));

        // redraw the data
        dataHist.plotOn(frame, MarkerStyle(kFullCircle), MarkerSize(0.8), DrawOption("ZP"), DataError(RooAbsData::SumW2));

        if (do_keys)
        {
            // rookeys

            lDataSet[ibin].Print();
            name.str("");
            name << "key_" << ibin;
            RooKeysPdf pdf_keys(name.str().c_str(), name.str().c_str(), lVar[ibin], lDataSet[ibin], RooKeysPdf::NoMirror, 2);

            RooPlot *xframe = lVar[ibin].frame(Bins(100), Title(Form("%s Wp_{T}=%.1f - %.1f GeV/c ", plabel, ptbins[ibin], ptbins[ibin + 2])));
            lDataSet[ibin].plotOn(xframe, DataError(RooAbsData::SumW2));
            TCanvas *c = new TCanvas("validatePDF", "validatePDF", 800, 600);
            c->cd();
            pdf_keys.plotOn(xframe, LineColor(kRed));
            xframe->Draw();
            c->SaveAs(Form("%s/plots/%s_%d_datasetW.pdf", outputDir, plabel, ibin));

            pdf_keys.plotOn(frame, LineColor(kRed));
            wksp->import(lDataSet[ibin]);
            wksp->import(pdf_keys);
            // wksp->Print();
        }

        int sizeParam = 0;
        if (string(plabel) == string("pfu1"))
            // sizeParam = 8; // 3 means + 3 sigma + 2 frac
            sizeParam = 4; // 1 means + 2 sigma + 1 frac
        if (string(plabel) == string("pfu2"))
            // sizeParam = 5; // 0 means + 3 sigma + 2 frac
            sizeParam = 3; // 0 means + 2 sigma + 1 frac

        TString nameRooHist = Form("h_%s", dataHist.GetName());
        TString nameRooCurve = Form("sig_%d_Norm[u_%d]", ibin, ibin);

        RooHist *hist = frame->getHist(nameRooHist.Data());
        RooCurve *fitCurve = frame->getCurve(nameRooCurve.Data());

        RooHist *hist_pull = hist->makePullHist(*fitCurve);
        hist_pull->SetTitle("");
        hist_pull->GetXaxis()->SetRangeUser(hv[ibin]->GetXaxis()->GetXmin(), hv[ibin]->GetXaxis()->GetXmax());
        hist_pull->GetXaxis()->SetTitle(xlabel);
        hist_pull->GetYaxis()->SetTitle("Pull");
        hist_pull->GetYaxis()->CenterTitle();
        hist_pull->GetYaxis()->SetRangeUser(-4.9, 4.9);
        hist_pull->GetYaxis()->SetNdivisions(511);
        hist_pull->SetMarkerColor(kAzure);
        hist_pull->SetLineColor(kAzure);
        hist_pull->SetFillColor(kAzure);
        //    hist_pull->GetYaxis()->SetTitleFont(42);
        //    hist_pull->GetXaxis()->SetTitleFont(42);
        hist_pull->GetYaxis()->SetTitleSize(0.045 * 3);
        hist_pull->GetYaxis()->SetTitleOffset(0.5);
        hist_pull->GetYaxis()->SetLabelOffset(0.014);
        hist_pull->GetYaxis()->SetLabelSize(0.040 * 3);
        hist_pull->GetYaxis()->SetLabelFont(42);
        hist_pull->GetXaxis()->SetTitleSize(0.045 * 3);
        hist_pull->GetXaxis()->SetTitleOffset(1.0);
        hist_pull->GetXaxis()->SetLabelOffset(0.014);
        hist_pull->GetXaxis()->SetLabelSize(0.040 * 3);
        hist_pull->GetXaxis()->SetLabelFont(42);

        double chi2Arr = frame->chiSquare(nameRooCurve.Data(), nameRooHist.Data(), sizeParam);
        double chi2ErrArr = 0;
        if (chi2Arr > 10)
        {
            chi2Arr = 10000;
            chi2ErrArr = 200;
        } // just a larger number so that is easy to notice on the plot
        //    cout << " chi2Arr=" << chi2Arr << " chi2ErrArr=" << chi2ErrArr << endl;

        if (sigOnly && do_5TeV)
            sprintf(lumi, "CMS Simumation                                                                 (5 TeV)");
        else if (sigOnly && !do_5TeV)
            sprintf(lumi, "CMS Simulation                                                                (13 TeV)");
        else if (!sigOnly && do_5TeV)
            sprintf(lumi, "CMS                                                              298.0 pb^{-1} (5 TeV)");
        else if (!sigOnly && !do_5TeV)
            sprintf(lumi, "CMS                                                             200.9 pb^{-1} (13 TeV)");
        sprintf(pname, "%sfit_%i", plabel, ibin);
        sprintf(ylabel, "Events / %.1f GeV", hv[ibin]->GetBinWidth(1));
        //    sprintf(binlabel,"%i < p_{T}(W) < %i",(Int_t)ptbins[ibin],(Int_t)ptbins[ibin+1]);
        sprintf(binlabel, "p_{T}(W) = %.1f - %.1f GeV ", ptbins[ibin], ptbins[ibin + 1]);

        sprintf(chi2text, "#chi^{2}/ndf = %.2f/%i", chi2Arr * ((Int_t)hv[ibin]->GetNbinsX() - sizeParam), (Int_t)hv[ibin]->GetNbinsX() - sizeParam);
        sprintf(hmeantext, "Mean = %.2f", hv[ibin]->GetMean());
        sprintf(hrmstext, "RMS = %.2f", hv[ibin]->GetRMS());

        if (etaBinCategory == 1)
            sprintf(binYlabel, "|y| < 0.5");
        if (etaBinCategory == 2)
            sprintf(binYlabel, "0.5 < |y| < 1");
        if (etaBinCategory == 3)
            sprintf(binYlabel, "|y| > 1");
        if (sigOnly)
        {
            // sprintf(nsigtext, "N_{evts} = %.1f", (Int_t)hv[ibin]->Integral());
            sprintf(nsigtext, "N_{sig} = %.1f #pm %.1f", nsig.getVal(), nsig.getError());
        }
        else
        {
            sprintf(nsigtext, "N_{sig} = %.1f #pm %.1f", nsig.getVal(), nsig.getError());
            sprintf(nbkgtext, "N_{bkg} = %.1f #pm %.1f", nbkg.getVal(), nbkg.getError());
        }
        sprintf(mean1text, "#mu_{1} = %.1f #pm %.1f", mean1Arr[ibin], mean1ErrArr[ibin]);
        sprintf(sig1text, "#sigma_{1} = %.1f #pm %.1f", sigma1Arr[ibin], sigma1ErrArr[ibin]);
        if (model >= 2)
        {
            sprintf(mean2text, "#mu_{2} = %.1f #pm %.1f", mean2Arr[ibin], mean2ErrArr[ibin]);
            // sprintf(sig0text,"#sigma = %.1f #pm %.1f",sigma0Arr[ibin],sigma0ErrArr[ibin]);
            sprintf(frac2text, "f_{2} = %.3f #pm %.3f", frac2Arr[ibin], frac2ErrArr[ibin]);
            sprintf(sig2text, "#sigma_{2} = %.1f #pm %.1f", sigma2Arr[ibin], sigma2ErrArr[ibin]);
        }
        if (model >= 3)
        {
            sprintf(mean3text, "#mu_{3} = %.1f #pm %.1f", mean3Arr[ibin], mean3ErrArr[ibin]);
            sprintf(sig3text, "#sigma_{3} = %.1f #pm %.1f", sigma3Arr[ibin], sigma3ErrArr[ibin]);
        }

        TCanvas *cLin = MakeCanvas("cLin", "cLin", 800, 800);
        cLin->Divide(1, 2, 0, 0);
        cLin->cd(1)->SetPad(0, 0.3, 1.0, 1.0);
        cLin->cd(1)->SetTopMargin(0.1);
        cLin->cd(1)->SetBottomMargin(0.01);
        cLin->cd(1)->SetLeftMargin(0.14);
        cLin->cd(1)->SetRightMargin(0.07);
        cLin->cd(1)->SetTickx(1);
        cLin->cd(1)->SetTicky(1);
        cLin->cd(2)->SetPad(0, 0, 1.0, 0.3);
        cLin->cd(2)->SetTopMargin(0.05);
        cLin->cd(2)->SetBottomMargin(0.40);
        cLin->cd(2)->SetLeftMargin(0.14);
        cLin->cd(2)->SetRightMargin(0.07);
        cLin->cd(2)->SetTickx(1);
        cLin->cd(2)->SetTicky(1);

        CPlot plot(pname, frame, "", xlabel, ylabel);
        plot.SetOutputDir(TString(outputDir) + "/plots");
        //    pad1->cd();
        plot.AddTextBox(lumi, 0.06, 0.92, 0.95, 0.97, 0, kBlack, -1);
        plot.AddTextBox(binlabel, 0.18, 0.80, 0.51, 0.85, 0, kBlack, -1);
        if (etaBinCategory != 0)
            plot.AddTextBox(binYlabel, 0.18, 0.71, 0.51, 0.66, 0, kBlack, -1);
        if (sigOnly)
            plot.AddTextBox(nsigtext, 0.18, 0.78, 0.51, 0.73, 0, kBlack, -1);
        //    else        plot.AddTextBox(0.18,0.78,0.51,0.68,0,kBlack,-1,2,nsigtext,nbkgtext);
        else
            plot.AddTextBox(nsigtext, 0.18, 0.78, 0.51, 0.73, 0, kBlack, -1); // this print the fraction now
        plot.AddTextBox(chi2text, 0.18, 0.71, 0.51, 0.66, 0, kBlack, -1);
        plot.AddTextBox(hmeantext, 0.18, 0.64, 0.51, 0.59, 0, kBlack, -1);
        plot.AddTextBox(hrmstext, 0.18, 0.57, 0.51, 0.52, 0, kBlack, -1);
        if (model == 1)
            plot.AddTextBox(0.70, 0.85, 0.95, 0.75, 0, kBlack, -1, 2, mean1text, sig1text);
        else if (model == 2)
            plot.AddTextBox(0.65, 0.83, 0.92, 0.60, 0, kBlack, -1, 4, mean1text, sig1text, sig2text, frac2text);
        // plot.AddTextBox(0.70, 0.85, 0.95, 0.65, 0, kBlack, -1, 5, mean1text, mean2text, sig1text, sig2text, frac2text);
        //    else if(model==3) plot.AddTextBox(0.70,0.90,0.95,0.65,0,kBlack,-1,7,mean1text,mean2text,mean3text,sig0text,sig1text,sig2text,sig3text);
        else if (model == 3)
            plot.AddTextBox(0.70, 0.85, 0.95, 0.60, 0, kBlack, -1, 6, mean1text, mean2text, mean3text, sig1text, sig2text, sig3text);
        // plot.Draw(cLin, kFALSE, "png", 1);
        plot.Draw(cLin, kFALSE, "pdf", 1);

        cLin->cd(2);
        hist_pull->Draw("A3 L ");
        TLine *lineZero = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0, hv[ibin]->GetXaxis()->GetXmax(), 0);
        lineZero->SetLineColor(kBlack);
        lineZero->Draw("same");
        TLine *lineZero1SigmaM = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0, hv[ibin]->GetXaxis()->GetXmax(), 0);
        lineZero1SigmaM->SetLineColor(11);
        lineZero1SigmaM->Draw("same");
        TLine *lineZero1SigmaP = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0, hv[ibin]->GetXaxis()->GetXmax(), 0);
        lineZero1SigmaP->SetLineColor(11);
        lineZero1SigmaP->Draw("same");

        // plot.Draw(cLin, kTRUE, "png", 1);
        plot.Draw(cLin, kTRUE, "pdf", 1);

        TCanvas *c1 = MakeCanvas("c1", "c1", 800, 800);
        c1->Divide(1, 2, 0, 0);
        c1->cd(1)->SetPad(0, 0.3, 1.0, 1.0);
        c1->cd(1)->SetTopMargin(0.1);
        c1->cd(1)->SetBottomMargin(0.01);
        c1->cd(1)->SetLeftMargin(0.14);
        c1->cd(1)->SetRightMargin(0.07);
        c1->cd(1)->SetTickx(1);
        c1->cd(1)->SetTicky(1);
        c1->cd(2)->SetPad(0, 0, 1.0, 0.3);
        c1->cd(2)->SetTopMargin(0.05);
        c1->cd(2)->SetBottomMargin(0.40);
        c1->cd(2)->SetLeftMargin(0.14);
        c1->cd(2)->SetRightMargin(0.07);
        c1->cd(2)->SetTickx(1);
        c1->cd(2)->SetTicky(1);

        sprintf(pname, "%sfitlog_%i", plabel, ibin);
        plot.SetYRange(0.1, 10 * hv[ibin]->GetMaximum());
        plot.SetName(pname);
        plot.SetLogy();
        // plot.Draw(c1, kFALSE, "png", 1);
        plot.Draw(c1, kFALSE, "pdf", 1);

        c1->cd(2);
        hist_pull->SetTitle("");
        hist_pull->Draw("A3 L ");
        lineZero->SetLineColor(kBlack);
        lineZero->Draw("same");
        lineZero1SigmaM->SetLineColor(11);
        lineZero1SigmaM->Draw("same");
        lineZero1SigmaP->SetLineColor(11);
        lineZero1SigmaP->Draw("same");
        // plot.Draw(c1, kTRUE, "png", 1);
        plot.Draw(c1, kTRUE, "pdf", 1);

        // reset color canvas
        c->SetFillColor(kWhite);
    }
}
