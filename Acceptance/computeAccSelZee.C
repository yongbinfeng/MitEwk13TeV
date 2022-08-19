//================================================================================================
//
// Compute Z->ee acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties no yet implemented
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TF1.h"
#include "TGraph.h"
#include "TLorentzVector.h" // 4-vector class
#include <TBenchmark.h> // class to track macro running statistics
#include <TClonesArray.h> // ROOT array class
#include <TFile.h> // file handle class
#include <TH1D.h> // histogram class
#include <TROOT.h> // access to gROOT, entry point to ROOT system
#include <TRandom3.h>
#include <TSystem.h> // interface to OS
#include <TTree.h> // class to access ntuples
#include <fstream> // functions for file I/O
#include <iomanip> // functions to format standard I/O
#include <iostream> // standard I/O
#include <sstream> // class for parsing strings
#include <string> // C++ string class
#include <vector> // STL vector class

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "MitEwk13TeV/EleScale/EnergyScaleCorrection.h" //EGMSmear
#include "MitEwk13TeV/Utils/LeptonCorr.hh" // Scale and resolution corrections

#include "MitEwk13TeV/Utils/CSample.hh" // helper class to handle samples
#include "MitEwk13TeV/Utils/ConfParse.hh" // input conf file parser
#include "MitEwk13TeV/Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "MitEwk13TeV/Utils/MyTools.hh" // various helper functions
#include "MitEwk13TeV/Utils/PrefiringEfficiency.cc" // prefiring efficiency functions

// helper class to handle efficiency tables
#include "MitEwk13TeV/Utils/CEffUser1D.hh"
#include "MitEwk13TeV/Utils/CEffUser2D.hh"
#include "MitEwk13TeV/Utils/AppEffSF.cc"
#endif

//=== MAIN MACRO =================================================================================================
void computeAccSelZee(const TString conf, // input file
    const TString inputDir,
    const TString outputDir, // output directory
    const TString outputName,
    const TString sqrts,
    const Int_t doPU,
    const Int_t doScaleCorr,
    const Int_t sigma,
    const TString SysFileGSFSel = "SysUnc_GSFSelEff.root",
    const bool is13TeV = 1)
{
    gBenchmark->Start("computeAccSelZee");
    const int gainSeed = 12;

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================

    const Double_t MASS_LOW = 60;
    const Double_t MASS_HIGH = 120;
    const Double_t PT_CUT = 25;
    // const Double_t ETA_CUT    = 2.5;
    const Double_t ETA_CUT = 2.4; // to match muons
    // const Double_t ETA_CUT    = 1.444; // to match muons
    const Double_t ELE_MASS = 0.000511;

    const Double_t ETA_BARREL = 1.4442;
    const Double_t ETA_ENDCAP = 1.566;
    // const Double_t ETA_BARREL = 10.;
    // const Double_t ETA_ENDCAP = 10.;

    const Int_t BOSON_ID = 23;
    const Int_t LEPTON_ID = 11;

    const int muEtaNB = 12;
    const float muEtaRange[muEtaNB + 1] = { -2.4, -2.0, -1.566, -1.4442, -1.0, -0.5, 0, 0.5, 1.0, 1.4442, 1.566, 2.0, 2.4 };
    const int muPtNB = 4;
    const float muPtRange[muPtNB + 1] = { 25, 30, 35, 40, 100000};

    const int NBptHLT = 12;
    const float ptrangeHLT[NBptHLT + 1] = { 25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 10000 };

    AppEffSF effs(inputDir);
    effs.loadHLT("EleHLTEff_aMCxPythia", "Positive", "Negative");
    effs.loadSel("EleGSFSelEff_aMCxPythia", "Combined", "Combined");
    effs.loadUncSel(SysFileGSFSel);

    //const TString corrFiles = "../EleScale/Run2017_LowPU_v2";
    const TString corrFiles = "/uscms_data/d3/yfeng/WpT/CMSSW_9_4_19/src/MitEwk13TeV/EleScale/Run2017_LowPU_v2";
    EnergyScaleCorrection eleCorr(corrFiles.Data(), EnergyScaleCorrection::ECALELF); // eleCorr.doScale= true; eleCorr.doSmearings =true;

    // load trigger menu
    const TString prefireEcalFileName = "/uscms/home/yfeng/nobackup/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Utils/All2017Gand2017HPrefiringMaps.root";
    const TString prefireMuonFileName = "/uscms/home/yfeng/nobackup/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Utils/L1MuonPrefiringParametriations.root";
    PrefiringEfficiency pfire(prefireEcalFileName.Data(), (is13TeV ? "2017H" : "2017G"), prefireMuonFileName.Data());

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================
    vector<TString> snamev; // sample name (for output files)
    vector<CSample*> samplev; // data/MC samples

    //
    // parse .conf file
    //
    confParse(conf, snamev, samplev);

    // Create output directory
    gSystem->mkdir(outputDir, kTRUE);

    TH2D* h = 0;

    TH2D* hHLTErr_pos = new TH2D("hHLTErr_pos", "", muEtaNB, muEtaRange, NBptHLT, ptrangeHLT);
    TH2D* hHLTErr_neg = new TH2D("hHLTErr_neg", "", muEtaNB, muEtaRange, NBptHLT, ptrangeHLT);

    TH2D* hGsfSelErr_pos = new TH2D("hGsfSelErr_pos", "", muEtaNB, muEtaRange, muPtNB, muPtRange);
    TH2D* hGsfSelErr_neg = new TH2D("hGsfSelErr_neg", "", muEtaNB, muEtaRange, muPtNB, muPtRange);

    // Data structures to store info from TTrees
    baconhep::TEventInfo* info = new baconhep::TEventInfo();
    baconhep::TGenEventInfo* gen = new baconhep::TGenEventInfo();
    TClonesArray* genPartArr = new TClonesArray("baconhep::TGenParticle");
    TClonesArray* electronArr = new TClonesArray("baconhep::TElectron");
    TClonesArray* vertexArr = new TClonesArray("baconhep::TVertex");
    TClonesArray* muonArr = new TClonesArray("baconhep::TMuon");
    TClonesArray* scArr = new TClonesArray("baconhep::TPhoton");
    TClonesArray* jetArr = new TClonesArray("baconhep::TJet");

    TFile* infile = 0;
    TTree* eventTree = 0;

    // Variables to store acceptances and uncertainties (per input file)
    Double_t nEvtsv = 0, nSelv = 0;
    Double_t nSelCorrv = 0, nSelCorrVarv = 0, nSelCorrVarvPos = 0, nSelCorrVarvNeg = 0;
    Double_t accv = 0, accCorrv = 0;
    Double_t accErrv = 0, accErrCorrv = 0, accErrCorrv_pos = 0, accErrCorrv_neg = 0;
    Double_t nSelCorrvFSR = 0, nSelCorrvMC = 0, nSelCorrvBkg = 0, nSelCorrvTag = 0;
    Double_t nSelCorrVarvFSR = 0, nSelCorrVarvMC = 0, nSelCorrVarvBkg = 0, nSelCorrVarvTag = 0;
    Double_t accCorrvFSR = 0, accCorrvMC = 0, accCorrvBkg = 0, accCorrvTag = 0;
    Double_t accErrCorrvFSR = 0, accErrCorrvMC = 0, accErrCorrvBkg = 0, accErrCorrvTag = 0;

    Double_t nSelPfire = 0, nSelPfireUp = 0, nSelPfireDown = 0;
    Double_t nSelPfireEcal = 0, nSelPfireEcalUp = 0, nSelPfireEcalDown = 0;
    Double_t nSelPfirePhoton = 0, nSelPfirePhotonUp = 0, nSelPfirePhotonDown = 0;
    Double_t nSelPfireJet = 0, nSelPfireJetUp = 0, nSelPfireJetDown = 0;
    Double_t nSelPfireMuon = 0, nSelPfireMuonUp = 0, nSelPfireMuonDown = 0;

    const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

    //
    // loop through files
    //
    CSample* samp = samplev[0];
    const UInt_t nfiles = samp->fnamev.size();

    for (UInt_t ifile = 0; ifile < nfiles; ifile++) {

        // Read input file and get the TTrees
        cout << "Processing " << samp->fnamev[ifile] << " ..." << endl;
        infile = TFile::Open(samp->fnamev[ifile]);
        assert(infile);

        eventTree = (TTree*)infile->Get("Events");
        assert(eventTree);
        eventTree->SetBranchAddress("Info", &info);
        TBranch* infoBr = eventTree->GetBranch("Info");
        eventTree->SetBranchAddress("GenEvtInfo", &gen);
        TBranch* genBr = eventTree->GetBranch("GenEvtInfo");
        eventTree->SetBranchAddress("GenParticle", &genPartArr);
        TBranch* genPartBr = eventTree->GetBranch("GenParticle");
        eventTree->SetBranchAddress("Electron", &electronArr);
        TBranch* electronBr = eventTree->GetBranch("Electron");
        eventTree->SetBranchAddress("PV", &vertexArr);
        TBranch* vertexBr = eventTree->GetBranch("PV");

        // jet and photon branches are 
        // for the prefire weights
        eventTree->SetBranchAddress("Photon", &scArr);
        TBranch* scBr = eventTree->GetBranch("Photon");
        eventTree->SetBranchAddress("AK4", &jetArr);
        TBranch* jetBr = eventTree->GetBranch("AK4");
        eventTree->SetBranchAddress("Muon", &muonArr);
        TBranch* muonBr = eventTree->GetBranch("Muon");


        //
        // loop over events
        //
        double frac = 0.01;
        if (!is13TeV)
            frac = 0.30;
        for (UInt_t ientry = 0; ientry < (uint)(frac * eventTree->GetEntries()); ientry++) {
            // for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry+=15) {
            if (ientry % 100000 == 0)
                cout << "Processing event " << ientry << ". " << (double)ientry / (double)eventTree->GetEntries() * 100 << " percent done with this file." << endl;
            genBr->GetEntry(ientry);
            genPartArr->Clear();
            genPartBr->GetEntry(ientry);
            infoBr->GetEntry(ientry);

            muonArr->Clear();
            muonBr->GetEntry(ientry);

            scArr->Clear();
            scBr->GetEntry(ientry);

            jetArr->Clear();
            jetBr->GetEntry(ientry);

            Int_t glepq1 = -99;
            Int_t glepq2 = -99;

            if (fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                continue;
            TLorentzVector* vec = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector* lep1 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector* lep2 = new TLorentzVector(0, 0, 0, 0);
            toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2, &glepq1, &glepq2, 1);
            //if (vec->M() < MASS_LOW || vec->M() > MASS_HIGH)
            //    continue;
            delete vec;
            delete lep1;
            delete lep2;

            vertexArr->Clear();
            vertexBr->GetEntry(ientry);
            double npv = vertexArr->GetEntries();
            Double_t weight = (gen->weight > 0) ? 1: -1;

            nEvtsv += weight;

            // trigger requirement
            if (!isEleTrigger(triggerMenu, info->triggerBits, kFALSE, is13TeV))
                continue;

            // good vertex requirement
            if (!(info->hasGoodPV))
                continue;

            electronArr->Clear();
            electronBr->GetEntry(ientry);

            for (Int_t i1 = 0; i1 < electronArr->GetEntriesFast(); i1++) {
                const baconhep::TElectron* ele1 = (baconhep::TElectron*)((*electronArr)[i1]);

                TLorentzVector vEle1(0, 0, 0, 0);
                vEle1.SetPtEtaPhiM(ele1->pt, ele1->eta, ele1->phi, ELE_MASS);

                if (doScaleCorr && (ele1->r9 < 1.)) {
                    // set up variable and apply smear correction to ele1
                    float ele1Smear = 0.;
                    float ele1Error = 0.;

                    float ele1AbsEta = fabs(vEle1.Eta());
                    float ele1Et = vEle1.E() / cosh(ele1AbsEta);

                    double tagEcalE = ele1->ecalEnergy;
                    double eTregress = tagEcalE / cosh(fabs(ele1->eta));
                    bool ele1isBarrel = ele1AbsEta < 1.4442;

                    float ele1R9Prime = ele1->r9; // r9 corrections MC only
                    double ele1Random = gRandom->Gaus(0, 1);

                    ele1Smear = eleCorr.smearingSigma(info->runNum, eTregress, ele1AbsEta, ele1R9Prime, gainSeed, 0., 0.);
                    float ele1SmearEP = eleCorr.smearingSigma(info->runNum, eTregress, ele1AbsEta, ele1R9Prime, gainSeed, 1., 0.);
                    float ele1SmearEM = eleCorr.smearingSigma(info->runNum, eTregress, ele1AbsEta, ele1R9Prime, gainSeed, -1., 0.);

                    if (sigma == 0)
                        (vEle1) *= 1. + ele1Smear * ele1Random;
                    else if (sigma == 1)
                        (vEle1) *= 1. + ele1SmearEP * ele1Random;
                    else if (sigma == -1)
                        (vEle1) *= 1. + ele1SmearEM * ele1Random;
                }

                if (vEle1.Pt() < PT_CUT)
                    continue; // lepton pT cut
                if (fabs(vEle1.Eta()) > ETA_CUT)
                    continue; // lepton |eta| cut
                if(fabs(vEle1.Eta())>=ETA_BARREL && fabs(vEle1.Eta())<=ETA_ENDCAP) 
                    continue;
                if (!passEleMediumID(ele1, vEle1, info->rhoIso))
                    continue; // lepton selection


                for (Int_t i2 = i1 + 1; i2 < electronArr->GetEntriesFast(); i2++) {
                    const baconhep::TElectron* ele2 = (baconhep::TElectron*)((*electronArr)[i2]);

                    TLorentzVector vEle2(0, 0, 0, 0);
                    vEle2.SetPtEtaPhiM(ele2->pt, ele2->eta, ele2->phi, ELE_MASS);

                    if (doScaleCorr && (ele2->r9 < 1.)) {
                        float ele2Smear = 0.;
                        float ele2Error = 0.;

                        float ele2AbsEta = fabs(vEle2.Eta());
                        float ele2Et = vEle2.E() / cosh(ele2AbsEta);
                        bool ele2isBarrel = ele2AbsEta < 1.4442;

                        float ele2R9Prime = ele2->r9; // r9 corrections MC only

                        double ele2e = ele2->ecalEnergy;
                        double eTregress2 = ele2e / cosh(fabs(ele2->eta));

                        double ele2Random = gRandom->Gaus(0, 1);

                        ele2Smear = eleCorr.smearingSigma(info->runNum, eTregress2, ele2AbsEta, ele2R9Prime, gainSeed, 0., 0.);
                        float ele2SmearEP = eleCorr.smearingSigma(info->runNum, eTregress2, ele2AbsEta, ele2R9Prime, gainSeed, 1., 0.);
                        float ele2SmearEM = eleCorr.smearingSigma(info->runNum, eTregress2, ele2AbsEta, ele2R9Prime, gainSeed, -1., 0.);

                        if (sigma == 0)
                            (vEle2) *= 1. + ele2Smear * ele2Random;
                        else if (sigma == 1)
                            (vEle2) *= 1. + ele2SmearEP * ele2Random;
                        else if (sigma == -1)
                            (vEle2) *= 1. + ele2SmearEM * ele2Random;
                    }

                    if (ele1->q == ele2->q)
                        continue;
                    if (vEle2.Pt() < PT_CUT)
                        continue; // lepton pT cut
                    if (fabs(vEle2.Eta()) > ETA_CUT)
                        continue; // lepton |eta| cut
                    if (!passEleMediumID(ele2, vEle2, info->rhoIso))
                        continue; // lepton selection
                    if(fabs(vEle2.Eta())>=ETA_BARREL && fabs(vEle2.Eta())<=ETA_ENDCAP) 
                        continue;

                    if (!isEleTriggerObj(triggerMenu, ele1->hltMatchBits, kFALSE, kFALSE, is13TeV) && !isEleTriggerObj(triggerMenu, ele2->hltMatchBits, kFALSE, kFALSE, is13TeV))
                        continue;

                    TLorentzVector vDilep = vEle1 + vEle2;
                    if ((vDilep.M() < MASS_LOW) || (vDilep.M() > MASS_HIGH))
                        continue;

                    /******** We have a Z candidate! HURRAY! ********/
                    Double_t effdata, effmc;
                    Double_t corr = 1;
                    Double_t effdataFSR, effdataMC, effdataBkg;

                    effdata = 1;
                    effmc = 1;
                    Double_t corrFSR = 1;
                    Double_t corrMC = 1;
                    Double_t corrBkg = 1;
                    Double_t corrTag = 1;
                    effdataFSR = 1;
                    effdataMC = 1;
                    effdataBkg = 1;
                    int q1 = ele1->q;
                    int q2 = ele2->q;

                    corr = effs.fullCorrections(&vEle1, q1, &vEle2, q2);
                    // corr = effs.dataOnly(&vEle1,q1,&vEle2,q2);
                    // dataeff = effs.dataOnly(&vEle1,q1,&vEle2,q2);
                    vector<double> uncs_gsf = effs.getUncSel(&vEle1, q1, &vEle2, q2);

                    corrFSR *= uncs_gsf[0] * effs.computeHLTSF(&vEle1, q1, &vEle2, q2); // alternate fsr model
                    corrMC *= uncs_gsf[1] * effs.computeHLTSF(&vEle1, q1, &vEle2, q2); // alternate mc gen model
                    corrBkg *= uncs_gsf[2] * effs.computeHLTSF(&vEle1, q1, &vEle2, q2); // alternate bkg model
                    corrTag *= uncs_gsf[3] * effs.computeHLTSF(&vEle1, q1, &vEle2, q2); // alternate bkg model
                    // corr *= effdata/effmc; // orig

                    double var = 0.;
                    // var += effs.statUncSta(&l1, q1) + effs.statUncSta(&l2, q2);
                    var += effs.statUncSel(&vEle1, q1, hGsfSelErr_pos, hGsfSelErr_neg, fabs(weight) * corr, true);
                    var += effs.statUncSel(&vEle2, q2, hGsfSelErr_pos, hGsfSelErr_neg, fabs(weight) * corr, true);
                    var += effs.statUncHLT(&vEle1, q1, hHLTErr_pos,    hHLTErr_neg,    fabs(weight) * corr);
                    var += effs.statUncHLT(&vEle2, q2, hHLTErr_pos,    hHLTErr_neg,    fabs(weight) * corr);
                    // cout << var1 << " " << var << endl;
                    // std::cout << "event " << info->evtNum << " weight " << corr << std::endl;
                    nSelv += weight;
                    nSelCorrvFSR += weight * corrFSR;
                    nSelCorrvMC += weight * corrMC;
                    nSelCorrvBkg += weight * corrBkg;
                    nSelCorrvTag += weight * corrTag;
                    nSelCorrv += weight * corr;
                    // std::cout << "corr " << corr << " corr FSR " << corrFSR << "  corr MC " << corrMC << "  corr Bkg " << corrBkg << std::endl;
                    nSelCorrVarvFSR += weight * weight * corrFSR * corrFSR;
                    nSelCorrVarvMC += weight * weight * corrMC * corrMC;
                    nSelCorrVarvBkg += weight * weight * corrBkg * corrBkg;
                    nSelCorrVarvTag += weight * weight * corrTag * corrTag;
                    // nSelCorrVarv+=weight*weight*corr*corr;

                    // prefire weights
                    float prefireEcal = 1, prefireEcalUp = 1, prefireEcalDown = 1;
                    float prefirePhoton = 1, prefirePhotonUp = 1, prefirePhotonDown = 1;
                    float prefireJet = 1, prefireJetUp = 1, prefireJetDown = 1;
                    float prefireMuon = 1, prefireMuonUp = 1, prefireMuonDown = 1, prefireMuonStatUp = 1, prefireMuonStatDown = 1, prefireMuonSystUp = 1, prefireMuonSystDown = 1;
                    pfire.setObjects(scArr, jetArr, muonArr);
                    pfire.computePhotonsOnly(prefirePhoton, prefirePhotonUp, prefirePhotonDown);
                    pfire.computeJetsOnly(prefireJet, prefireJetUp, prefireJetDown);
                    pfire.computeEcalsOnly(prefireEcal, prefireEcalUp, prefireEcalDown);
                    pfire.computeMuonsOnly(prefireMuon, prefireMuonUp, prefireMuonDown, prefireMuonStatUp, prefireMuonStatDown, prefireMuonSystUp, prefireMuonSystDown);
                    double prefireWeight = prefireEcal * prefireMuon;
                    double prefireUp = prefireEcalUp * prefireMuonUp;
                    double prefireDown = prefireEcalDown * prefireMuonDown;

                    nSelPfire         += weight * corr * prefireWeight;
                    nSelPfireUp       += weight * corr * prefireUp;
                    nSelPfireDown     += weight * corr * prefireDown;
                    // ecal-prefire-related
                    nSelPfireEcal     += weight * corr * prefireEcal;
                    nSelPfireEcalUp   += weight * corr * prefireEcalUp;
                    nSelPfireEcalDown += weight * corr * prefireEcalDown;
                    // muon-prefire-related
                    nSelPfireMuon     += weight * corr * prefireMuon;
                    nSelPfireMuonUp   += weight * corr * prefireMuonUp;
                    nSelPfireMuonDown += weight * corr * prefireMuonDown;
                    // photon-prefire-related
                    nSelPfirePhoton   += weight * corr * prefirePhoton;
                    nSelPfirePhotonUp += weight * corr * prefirePhotonUp;
                    nSelPfirePhotonDown += weight * corr * prefirePhotonDown;
                    // jet-prefire-related
                    nSelPfireJet      += weight * corr * prefireJet;
                    nSelPfireJetUp    += weight * corr * prefireJetUp;
                    nSelPfireJetDown  += weight * corr * prefireJetDown;

                }
            }
        }

        Double_t var = 0, var_pos = 0, var_neg = 0;
        for (Int_t iy = 0; iy <= hHLTErr_pos->GetNbinsY() + 1; iy++) {
            for (Int_t ix = 0; ix <= hHLTErr_pos->GetNbinsX() + 1; ix++) {
                Double_t err = hHLTErr_pos->GetBinContent(ix, iy);
                var += err * err;
                var_pos += err * err;

                err = hHLTErr_neg->GetBinContent(ix, iy);
                var += err * err;
                var_neg += err * err;
            }
        }
        for (Int_t iy = 0; iy <= hGsfSelErr_pos->GetNbinsY() + 1; iy++) {
            for (Int_t ix = 0; ix <= hGsfSelErr_pos->GetNbinsX() + 1; ix++) {
                Double_t err = hGsfSelErr_pos->GetBinContent(ix, iy);
                var += err * err;
                var_pos += err * err;

                err = hGsfSelErr_neg->GetBinContent(ix, iy);
                var += err * err;
                var_neg += err * err;
            }
        }
        std::cout << "blah  " << std::endl;
        nSelCorrVarvFSR += var;
        nSelCorrVarvMC += var;
        nSelCorrVarvBkg += var;
        nSelCorrVarv += var;
        nSelCorrVarvPos += var_pos;
        nSelCorrVarvNeg += var_neg;
        cout << var << endl;
        std::cout << "comput acceptances" << std::endl;

        // compute acceptances
        std::cout << nEvtsv << " " << nSelv << std::endl;
        accv = (nSelv / nEvtsv);

        accCorrvFSR = (nSelCorrvFSR / nEvtsv);
        accCorrvMC = (nSelCorrvMC / nEvtsv);
        accCorrvBkg = (nSelCorrvBkg / nEvtsv);
        accCorrvTag = (nSelCorrvTag / nEvtsv);
        accCorrv = (nSelCorrv / nEvtsv);

        accErrv = (sqrt(accv * (1. + accv) / nEvtsv));
        accErrCorrvFSR = (accCorrvFSR * sqrt((nSelCorrVarv) / (nSelCorrvFSR * nSelCorrvFSR) + 1. / nEvtsv));
        accErrCorrvMC = (accCorrvMC * sqrt((nSelCorrVarv) / (nSelCorrvMC * nSelCorrvMC) + 1. / nEvtsv));
        accErrCorrvBkg = (accCorrvBkg * sqrt((nSelCorrVarv) / (nSelCorrvBkg * nSelCorrvBkg) + 1. / nEvtsv));
        accErrCorrvTag = (accCorrvTag * sqrt((nSelCorrVarv) / (nSelCorrvTag * nSelCorrvTag) + 1. / nEvtsv));
        //accErrCorrv = (accCorrv * sqrt((nSelCorrVarv) / (nSelCorrv * nSelCorrv) + 1. / nEvtsv));
        //accErrCorrv_pos = (accCorrv * sqrt((nSelCorrVarvPos) / (nSelCorrv * nSelCorrv) + 1. / nEvtsv));
        //accErrCorrv_neg = (accCorrv * sqrt((nSelCorrVarvNeg) / (nSelCorrv * nSelCorrv) + 1. / nEvtsv));
        accErrCorrv = (accCorrv * sqrt(nSelCorrVarv) / nSelCorrv );
        accErrCorrv_pos = (accCorrv * sqrt(nSelCorrVarvPos) / nSelCorrv );
        accErrCorrv_neg = (accCorrv * sqrt(nSelCorrVarvNeg) / nSelCorrv );

        delete infile;
        infile = 0, eventTree = 0;
    }
    delete info;
    delete gen;
    delete electronArr;

    //--------------------------------------------------------------------------------------------------------------
    // Output
    //==============================================================================================================

    // Print full set for efficiency calculations
    char masterOutput[600];
    // just start printing....
    //for (uint ifile = 0; ifile < fnamev.size(); ++ifile) { // go through info per file
    for (uint ifile = 0; ifile < 1; ++ifile) { // go through info per file
        sprintf(masterOutput, "%s/%s.txt", outputDir.Data(), outputName.Data());
        ofstream txtfile;
        txtfile.open(masterOutput);
        txtfile << "acc " << accCorrv << endl;
        txtfile << "acc_stat " << accCorrv + accErrCorrv << endl;
        txtfile << "acc_stat pos " << accCorrv + accErrCorrv_pos << endl;
        txtfile << "acc_stat neg " << accCorrv + accErrCorrv_neg << endl;
        txtfile << "sel_fsr"
                << " " << accCorrvFSR << endl;
        txtfile << "sel_mc"
                << " " << accCorrvMC << endl;
        txtfile << "sel_bkg"
                << " " << accCorrvBkg << endl;
        txtfile << "sel_tagpt"
                << " " << accCorrvTag << endl;
        txtfile.close();
    }

    double acc_sys = 0.;
    int ifile = 0;
    acc_sys += TMath::Power(accErrCorrv_pos / accCorrv, 2.0);
    acc_sys += TMath::Power(accErrCorrv_neg / accCorrv, 2.0);
    acc_sys += TMath::Power(accCorrvFSR / accCorrv - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvMC / accCorrv - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvBkg / accCorrv - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvTag / accCorrv - 1.0, 2.0);
    acc_sys = TMath::Sqrt(acc_sys);

    cout << "*" << endl;
    cout << "* SUMMARY" << endl;
    cout << "*--------------------------------------------------" << endl;
    cout << " Z -> e e" << endl;
    cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
    cout << "  pT > " << PT_CUT << endl;
    cout << "  |eta| < " << ETA_CUT << endl;
    cout << endl;

    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        cout << "    *** Acceptance ***" << endl;
        cout << "          nominal: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv << endl;
        cout << "     SF corrected: " << accCorrv << " +/- " << accErrCorrv << endl;
        cout << "     SF corrected pos: " << accCorrv << " +/- " << accErrCorrv_pos << endl;
        cout << "     SF corrected neg: " << accCorrv << " +/- " << accErrCorrv_neg << endl;
        cout << "  ==total efficiency==> " << setw(4) << accCorrv / accv << endl;
        cout << "          stat: " << 100 * accErrCorrv / accCorrv << endl;
        cout << "          stat pos: " << 100 * accErrCorrv_pos / accCorrv << endl;
        cout << "          stat neg: " << 100 * accErrCorrv_neg / accCorrv << endl;
        cout << "          FSR unc: " << 100 * (accCorrvFSR / accCorrv - 1.0) << " +/- " << accErrCorrvFSR << endl;
        cout << "           MC unc: " << 100 * (accCorrvMC / accCorrv - 1.0) << " +/- " << accErrCorrvMC << endl;
        cout << "          Bkg unc: " << 100 * (accCorrvBkg / accCorrv - 1.0) << " +/- " << accErrCorrvBkg << endl;
        cout << "          Tag unc: " << 100 * (accCorrvTag / accCorrv - 1.0) << " +/- " << accErrCorrvTag << endl;
        cout << "          Total: " << acc_sys << endl;
        cout << endl;

        cout << endl << endl;
        cout << " Prefire Correction " << nSelPfire / nSelCorrv << " + " << fabs(nSelPfireUp - nSelPfire) / nSelPfire << " - " << fabs(nSelPfireDown - nSelPfire) / nSelPfire << endl;
        cout << " Prefire ECAL Correction " << nSelPfireEcal / nSelCorrv << " + " << fabs(nSelPfireEcalUp - nSelPfireEcal) / nSelPfireEcal << " - " << fabs(nSelPfireEcalDown - nSelPfireEcal) / nSelPfireEcal << endl;
        cout << " Prefire Muon Correction " << nSelPfireMuon / nSelCorrv << " + " << fabs(nSelPfireMuonUp - nSelPfireMuon) / nSelPfireMuon << " - " << fabs(nSelPfireMuonDown - nSelPfireMuon) / nSelPfireMuon << endl;
        cout << " Prefire Photon Correction " << nSelPfirePhoton / nSelCorrv << " + " << fabs(nSelPfirePhotonUp - nSelPfirePhoton) / nSelPfirePhoton << " - " << fabs(nSelPfirePhotonDown - nSelPfirePhoton) / nSelPfirePhoton << endl;
        cout << " Prefire Jet Correction " << nSelPfireJet / nSelCorrv << " + " << fabs(nSelPfireJetUp - nSelPfireJet) / nSelPfireJet << " - " << fabs(nSelPfireJetDown - nSelPfireJet) / nSelPfireJet << endl;


    }

    char txtfname[500];
    sprintf(txtfname, "%s/binned.txt", outputDir.Data());
    ofstream txtfile;
    txtfile.open(txtfname);
    txtfile << "*" << endl;
    txtfile << "* SUMMARY" << endl;
    txtfile << "*--------------------------------------------------" << endl;
    txtfile << " Z -> e e" << endl;
    txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
    txtfile << "  pT > " << PT_CUT << endl;
    txtfile << "  |eta| < " << ETA_CUT << endl;
    txtfile << endl;

    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        txtfile << "    *** Acceptance ***" << endl;
        txtfile << "          nominal: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv << endl;
        txtfile << "  ==total efficiency==> " << setw(4) << accCorrv / accv << endl;
        txtfile << "     SF corrected: " << accCorrv << " +/- " << accErrCorrv << endl;
        txtfile << "     SF corrected Pos: " << accCorrv << " +/- " << accErrCorrv_pos << endl;
        txtfile << "     SF corrected Neg: " << accCorrv << " +/- " << accErrCorrv_neg << endl;
        txtfile << "          FSR unc: " << accCorrvFSR << " +/- " << accErrCorrvFSR << endl;
        txtfile << "           MC unc: " << accCorrvMC << " +/- " << accErrCorrvMC << endl;
        txtfile << "          Bkg unc: " << accCorrvBkg << " +/- " << accErrCorrvBkg << endl;
        txtfile << "          Bkg unc: " << accCorrvTag << " +/- " << accErrCorrvTag << endl;
        txtfile << endl;
    }
    txtfile.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sel_nums_only.txt", outputDir.Data());
    ofstream txtfile2;
    txtfile2.open(txtfname);

    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        txtfile2 << accCorrv << " " << accErrCorrv << endl;
        txtfile2 << accCorrvFSR << endl;
        txtfile2 << accCorrvMC << endl;
        txtfile2 << accCorrvBkg << endl;
        txtfile2 << accCorrvTag << endl;
        // txtfile << accCorrvFSR  << ", " << accCorrvMC << ", " << accCorrvBkg << ", " << accCorrvTag << endl;

        txtfile2 << endl;
    }
    txtfile2.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/gsf_unc.txt", outputDir.Data());
    ofstream txtfile3;
    txtfile3.open(txtfname);

    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        txtfile3 << accCorrv << endl;
        txtfile3 << accCorrvFSR << endl;
        txtfile3 << accCorrvMC << endl;
        txtfile3 << accCorrvBkg << endl;
        txtfile3 << accCorrvTag << endl;
        // txtfile << accCorrvFSR  << ", " << accCorrvMC << ", " << accCorrvBkg << ", " << accCorrvTag << endl;

        txtfile3 << endl;
    }
    txtfile3.close();

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;

    gBenchmark->Show("computeAccSelZee");
}
