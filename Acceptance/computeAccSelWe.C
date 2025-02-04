//================================================================================================
//
// Compute W->enu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties need to be checked
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TF1.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include <TBenchmark.h> // class to track macro running statistics
#include <TClonesArray.h> // ROOT array class
#include <TFile.h> // file handle class
#include <TH1D.h> // histogram class
#include <TMath.h> // ROOT math library
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
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "MitEwk13TeV/Utils/LeptonCorr.hh" // Scale and resolution corrections
#include "MitEwk13TeV/Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "MitEwk13TeV/Utils/MyTools.hh"
#include "MitEwk13TeV/EleScale/EnergyScaleCorrection.h" //EGMSmear
#include "MitEwk13TeV/Utils/PrefiringEfficiency.cc" // prefiring efficiency functions

#include "MitEwk13TeV/Utils/ConfParse.hh" // input conf file parser
#include "MitEwk13TeV/Utils/CSample.hh" // helper class to handle samples
// helper class to handle efficiency tables
#include "MitEwk13TeV/Utils/AppEffSF.cc"
#endif

//=== MAIN MACRO =================================================================================================

void computeAccSelWe(const TString conf, // input file
    const TString inputDir, // efficiency main directory
    const TString outputDir, // output directory
    const TString outputName,
    const Int_t charge, // 0 = inclusive, +1 = W+, -1 = W-
    const Int_t doPU,
    const Int_t doScaleCorr,
    const Int_t sigma,
    const TString SysFileGSFSel = "SysUnc_GSFSelEff.root",
    const bool is13TeV = 1)
{
    gBenchmark->Start("computeAccSelWe");
    const int gainSeed = 12;
    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================

    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;
    const Double_t ELE_MASS = 0.000511;

    const Double_t ETA_BARREL = 1.4442;
    const Double_t ETA_ENDCAP = 1.566;
    // const Double_t ETA_BARREL = 10.;
    // const Double_t ETA_ENDCAP = 10.;

    const Double_t VETO_PT = 10;
    const Double_t VETO_ETA = 2.4;

    const Int_t BOSON_ID = 24;
    const Int_t LEPTON_ID = 11;

    const int muEtaNB = 12;
    const float muEtaRange[muEtaNB + 1] = { -2.4, -2.0, -1.566, -1.4442, -1.0, -0.5, 0, 0.5, 1.0, 1.4442, 1.566, 2.0, 2.4 };
    const int muPtNB = 4;
    const float muPtRange[muPtNB + 1] = { 25, 30, 35, 40, 1000000};

    const int NBptHLT = 12;
    const float ptrangeHLT[NBptHLT + 1] = { 25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 10000 };

    AppEffSF effs(inputDir);
    effs.loadHLT("EleHLTEff_aMCxPythia", "Positive", "Negative");
    effs.loadSel("EleGSFSelEff_aMCxPythia", "Combined", "Combined");
    effs.loadUncSel(SysFileGSFSel);

    const TString corrFiles = "/uscms/home/yfeng/nobackup/WpT/CMSSW_9_4_19/src/MitEwk13TeV/EleScale/Run2017_LowPU_v2";

    EnergyScaleCorrection eleCorr(corrFiles.Data(), EnergyScaleCorrection::ECALELF);

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

    TH2D* hHLTErr_pos = new TH2D("hHLTErr_pos", "", muEtaNB, muEtaRange, NBptHLT, ptrangeHLT);
    TH2D* hHLTErr_neg = new TH2D("hHLTErr_neg", "", muEtaNB, muEtaRange, NBptHLT, ptrangeHLT);

    TH2D* hGsfSelErr_pos = new TH2D("hGsfSelErr_pos", "", muEtaNB, muEtaRange, muPtNB, muPtRange);
    TH2D* hGsfSelErr_neg = new TH2D("hGsfSelErr_neg", "", muEtaNB, muEtaRange, muPtNB, muPtRange);

    // Data structures to store info from TTrees
    baconhep::TEventInfo* info = new baconhep::TEventInfo();
    baconhep::TGenEventInfo* gen = new baconhep::TGenEventInfo();
    TClonesArray* electronArr = new TClonesArray("baconhep::TElectron");
    TClonesArray* genPartArr = new TClonesArray("baconhep::TGenParticle");
    TClonesArray* vertexArr = new TClonesArray("baconhep::TVertex");
    TClonesArray* scArr = new TClonesArray("baconhep::TPhoton");
    TClonesArray* jetArr = new TClonesArray("baconhep::TJet");
    TClonesArray* muonArr = new TClonesArray("baconhep::TMuon");

    TFile* infile = 0;
    TTree* eventTree = 0;

    // Variables to store acceptances and uncertainties
    Double_t nEvtsv = 0, nSelv = 0, nSelBv = 0, nSelEv = 0;
    Double_t accv = 0, accBv = 0, accEv = 0;
    Double_t accErrv = 0, accErrBv = 0, accErrEv = 0;
    Double_t nSelCorrv = 0, nSelBCorrv = 0, nSelECorrv = 0;
    Double_t nSelCorrVarv = 0, nSelBCorrVarv = 0, nSelECorrVarv = 0, nSelCorrVarvPos = 0, nSelCorrVarvNeg = 0;
    Double_t accCorrv = 0, accBCorrv = 0, accECorrv = 0;
    Double_t accErrCorrv = 0, accErrBCorrv = 0, accErrECorrv = 0, accErrCorrv_pos = 0, accErrCorrv_neg = 0;
    Double_t nSelCorrvFSR = 0, nSelCorrvMC = 0, nSelCorrvBkg = 0, nSelCorrvTag = 0;
    Double_t accCorrvFSR = 0, accCorrvMC = 0, accCorrvBkg = 0, accCorrvTag = 0;
    Double_t accErrCorrvFSR = 0, accErrCorrvMC = 0, accErrCorrvBkg = 0, accErrCorrvTag = 0;

    Double_t nSelPfire = 0, nSelPfireUp = 0, nSelPfireDown = 0;
    Double_t nSelPfireEcal = 0, nSelPfireEcalUp = 0, nSelPfireEcalDown = 0;
    Double_t nSelPfirePhoton = 0, nSelPfirePhotonUp = 0, nSelPfirePhotonDown = 0;
    Double_t nSelPfireJet = 0, nSelPfireJetUp = 0, nSelPfireJetDown = 0;
    Double_t nSelPfireMuon = 0, nSelPfireMuonUp = 0, nSelPfireMuonDown = 0;

    const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

    //
    // loop through samples
    //
    for (UInt_t isamp = 0; isamp < samplev.size(); isamp++) {

        CSample* samp = samplev[isamp];
        const UInt_t nfiles = samp->fnamev.size();

        // loop through files
        for (UInt_t ifile = 0; ifile < nfiles; ifile++) {

            // Read input file and get the TTrees
            cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... ";
            cout.flush();
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
            double frac = 0.05; // fraction of events to be used for calculation
            if (isamp == 0) frac = 0.05;
            if (isamp == 1) frac = 0.10;
            if (isamp == 2) frac = 0.20;
            if (!is13TeV) frac = 0.30;
            double nWgtSum = 0., nAbsSum = 0; // total number of events after reweighting

            // loop over the events first, to get the positive and negative frations of events,
            // used for scaling later
            std::cout << "Process events quickly first time for counting" << std::endl;
            for (UInt_t ientry = 0; ientry < (uint)(frac * eventTree->GetEntries()); ientry++) {
                if (ientry % 100000 == 0)
                    cout << "Processing event " << ientry << ". " << (double)ientry / (double)eventTree->GetEntries() * 100 << " percent done with this file." << endl;
                genBr->GetEntry(ientry);
                nAbsSum += 1.;
                nWgtSum += (gen->weight > 0) ? 1 : -1;
            }
            std::cout << "Finished first loop. Total events " << nAbsSum << " after negative weight subtraction " << nWgtSum / nAbsSum << std::endl;

            for (UInt_t ientry = 0; ientry < (uint)(frac * eventTree->GetEntries()); ientry++) {
                if (ientry % 100000 == 0)
                    cout << "Processing event " << ientry << ". " << (double)ientry / (double)eventTree->GetEntries() * 100 << " percent done with this file." << endl;
                infoBr->GetEntry(ientry);
                genBr->GetEntry(ientry);
                genPartArr->Clear();
                genPartBr->GetEntry(ientry);

                if (charge == -1 && toolbox::flavor(genPartArr, -BOSON_ID) != LEPTON_ID)
                    continue;
                if (charge == 1 && toolbox::flavor(genPartArr, BOSON_ID) != -LEPTON_ID)
                    continue;
                if (charge == 0 && fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                    continue;
                /*TLorentzVector *vec=new TLorentzVector(0,0,0,0);
                TLorentzVector *lep1=new TLorentzVector(0,0,0,0);
                TLorentzVector *lep2=new TLorentzVector(0,0,0,0);
                toolbox::fillGen(genPartArr, BOSON_ID, vec, lep1, lep2,1);*/
                Int_t glepq1 = -99;
                Int_t glepq2 = -99;
                TLorentzVector* gvec = new TLorentzVector(0, 0, 0, 0);
                TLorentzVector* glep1 = new TLorentzVector(0, 0, 0, 0);
                TLorentzVector* glep2 = new TLorentzVector(0, 0, 0, 0);
                toolbox::fillGen(genPartArr, BOSON_ID, gvec, glep1, glep2, &glepq1, &glepq2, 1);

                // TLorentzVector tvec=*glep1+*glep2;
                // TLorentzVector* genV=new TLorentzVector(0,0,0,0);
                // genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
                // genVPt   = tvec.Pt();
                // genVPhi  = tvec.Phi();
                // genVy    = tvec.Rapidity();
                // double genVMass = tvec.M();
                //double mtgen = sqrt(2.0 * (glep1->Pt()) * (glep2->Pt()) * (1.0 - cos(toolbox::deltaPhi(glep1->Phi(), glep2->Phi()))));
                // if(mtgen < 40) continue;
                // if(mtgen < 40 || mtgen > 140) continue;
                // cout << "mass " << genVMass <<  "   mt " << mt << endl;

                vertexArr->Clear();
                vertexBr->GetEntry(ientry);
                double npv = vertexArr->GetEntries();
                Double_t weight = (gen->weight > 0 ? 1 : -1) * samp->xsecv[ifile] / nWgtSum;

                nEvtsv += weight;

                if (!isEleTrigger(triggerMenu, info->triggerBits, kFALSE, is13TeV))
                    continue;

                // good vertex requirement
                if (!(info->hasGoodPV))
                    continue;

                electronArr->Clear();
                electronBr->GetEntry(ientry);

                muonArr->Clear();
                muonBr->GetEntry(ientry);

                scArr->Clear();
                scBr->GetEntry(ientry);

                jetArr->Clear();
                jetBr->GetEntry(ientry);

                Int_t nLooseLep = 0;
                const baconhep::TElectron* goodEle = 0;
                TLorentzVector vEle(0, 0, 0, 0);
                TLorentzVector vElefinal(0, 0, 0, 0);
                Bool_t passSel = kFALSE;
                int q = 0;
                for (Int_t i = 0; i < electronArr->GetEntriesFast(); i++) {
                    const baconhep::TElectron* ele = (baconhep::TElectron*)((*electronArr)[i]);
                    vEle.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, ELE_MASS);
                    // check ECAL gap
                    if (doScaleCorr && (ele->r9 < 1.)) {
                        float eleSmear = 0.;

                        float eleAbsEta = fabs(vEle.Eta());
                        float eleEt = vEle.E() / cosh(eleAbsEta);
                        bool eleisBarrel = eleAbsEta < 1.4442;

                        float eleR9Prime = ele->r9; // r9 corrections MC only, none after 2016

                        double eleRamdom = gRandom->Gaus(0, 1);

                        eleSmear = eleCorr.smearingSigma(info->runNum, eleEt, eleAbsEta, eleR9Prime, gainSeed, 0., 0.);
                        float eleSmearEP = eleCorr.smearingSigma(info->runNum, eleEt, eleAbsEta, eleR9Prime, gainSeed, 1., 0.);
                        float eleSmearEM = eleCorr.smearingSigma(info->runNum, eleEt, eleAbsEta, eleR9Prime, gainSeed, -1., 0.);

                        if (sigma == 0)
                            (vEle) *= 1. + eleSmear * eleRamdom;
                        else if (sigma == 1)
                            (vEle) *= 1. + eleSmearEP * eleRamdom;
                        else if (sigma == -1)
                            (vEle) *= 1. + eleSmearEM * eleRamdom;
                    }

                    if (fabs(vEle.Eta()) > VETO_ETA)
                        continue;
                    if (vEle.Pt() < VETO_PT)
                        continue;
                    if (passEleLooseID(ele, vEle, info->rhoIso))
                        nLooseLep++;

                    if (nLooseLep > 1) { // extra lepton veto
                        passSel = kFALSE;
                        break;
                    }

                    if (vEle.Pt() < PT_CUT)
                        continue; // lepton pT cut
                    if (fabs(vEle.Eta()) > ETA_CUT)
                        continue; // lepton |eta| cut
                    if(fabs(vEle.Eta())>=ETA_BARREL && fabs(vEle.Eta())<=ETA_ENDCAP)
                        continue;
                    if (!passEleMediumID(ele, vEle, info->rhoIso))
                        continue; // lepton selection

                    if (!isEleTriggerObj(triggerMenu, ele->hltMatchBits, kFALSE, kFALSE, is13TeV))
                        continue;
                    q = ele->q;
                    if (charge != 0 && ele->q != charge)
                        continue; // check charge (if necessary)

                    double mtreco = sqrt(2.0 * (ele->pt) * (info->pfMETC) * (1.0 - cos(toolbox::deltaPhi(ele->phi, info->pfMETCphi))));

                    if (mtreco < 20)
                        continue;
                    // if(mtreco > 40 && mtreco < 140) continue;

                    passSel = kTRUE;
                    goodEle = ele;
                    vElefinal = vEle;
                }

                if (passSel) {

                    /******** We have a W candidate! HURRAY! ********/
                    Bool_t isBarrel = (fabs(vElefinal.Eta()) < ETA_BARREL) ? kTRUE : kFALSE;

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

                    corr = effs.fullCorrections(&vElefinal, q);
                    // corr = effs.dataOnly(&vElefinal,q);
                    vector<double> uncs_gsf = effs.getUncSel(&vElefinal, q);

                    corrFSR *= uncs_gsf[0] * effs.computeHLTSF(&vElefinal, q); // alternate fsr model
                    corrMC *= uncs_gsf[1] * effs.computeHLTSF(&vElefinal, q); // alternate mc gen model
                    corrBkg *= uncs_gsf[2] * effs.computeHLTSF(&vElefinal, q); // alternate bkg model
                    corrTag *= uncs_gsf[3] * effs.computeHLTSF(&vElefinal, q); // alternate bkg model

                    double var = 0.;
                    var += effs.statUncSel(&vElefinal, q, hGsfSelErr_pos, hGsfSelErr_neg, fabs(weight) * corr, true);
                    var += effs.statUncHLT(&vElefinal, q, hHLTErr_pos,    hHLTErr_neg,    fabs(weight) * corr);

                    nSelv += weight;
                    nSelCorrv += weight * corr;
                    // nSelCorrVarv[ifile]+=weight*weight*corr*corr;
                    nSelCorrvFSR += weight * corrFSR;
                    nSelCorrvMC += weight * corrMC;
                    nSelCorrvBkg += weight * corrBkg;
                    nSelCorrvTag += weight * corrTag;

                    if (isBarrel) {
                        nSelBv += weight;
                        nSelBCorrv += weight * corr;
                        nSelBCorrVarv += weight * weight * corr * corr;

                    } else {
                        nSelEv += weight;
                        nSelECorrv += weight * corr;
                        nSelECorrVarv += weight * weight * corr * corr;
                    }

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

            //delete infile;
            //infile = 0, eventTree = 0;

            std::cout << " nSelCorrVarv " << nSelCorrVarv << " nSelCorrVarvPos " << nSelCorrVarvPos << " nSelCorrVarvNeg " << nSelCorrVarvNeg << std::endl;
        }
    }


    Double_t var = 0, varB = 0, varE = 0, var_pos = 0, var_neg = 0;
    for (Int_t iy = 0; iy <= hHLTErr_pos->GetNbinsY() + 1; iy++) {
        for (Int_t ix = 0; ix <= hHLTErr_pos->GetNbinsX() + 1; ix++) {
            Double_t err = hHLTErr_pos->GetBinContent(ix, iy);
            var += err * err;
            var_pos += err * err;
            // cout << "hlt pos " << err*err << endl;
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
            // cout << "gsf pos " << err*err << endl;
        }
    }

    nSelCorrVarv += var;
    nSelCorrVarvPos += var_pos;
    nSelCorrVarvNeg += var_neg;
    nSelBCorrVarv += varB;
    nSelECorrVarv += varE;

    // compute acceptances
    accv = (nSelv / nEvtsv);
    accErrv = (sqrt(accv * (1. + accv) / nEvtsv));
    accBv = (nSelBv / nEvtsv);
    accErrBv = (sqrt(accBv * (1. + accBv) / nEvtsv));
    accEv = (nSelEv / nEvtsv);
    accErrEv = (sqrt(accEv * (1. + accEv) / nEvtsv));

    accCorrv = (nSelCorrv / nEvtsv);
    //accErrCorrv = (accCorrv * sqrt(nSelCorrVarv / nSelCorrv / nSelCorrv + 1. / nEvtsv));
    //accErrCorrv_pos = (accCorrv * sqrt(nSelCorrVarvPos / nSelCorrv / nSelCorrv + 1. / nEvtsv));
    //accErrCorrv_neg = (accCorrv * sqrt(nSelCorrVarvNeg / nSelCorrv / nSelCorrv + 1. / nEvtsv));
    accErrCorrv = (accCorrv * sqrt(nSelCorrVarv) / nSelCorrv );
    accErrCorrv_pos = (accCorrv * sqrt(nSelCorrVarvPos) / nSelCorrv );
    accErrCorrv_neg = (accCorrv * sqrt(nSelCorrVarvNeg) / nSelCorrv );
    accBCorrv = (nSelBCorrv / nEvtsv);
    accErrBCorrv = (accBCorrv * sqrt(nSelBCorrVarv / nSelBCorrv / nSelBCorrv + 1. / nEvtsv));
    accECorrv = (nSelECorrv / nEvtsv);
    accErrECorrv = (accECorrv * sqrt(nSelECorrVarv / nSelECorrv / nSelECorrv + 1. / nEvtsv));

    accCorrvFSR = (nSelCorrvFSR / nEvtsv);
    accCorrvMC = (nSelCorrvMC / nEvtsv);
    accCorrvBkg = (nSelCorrvBkg / nEvtsv);
    accCorrvTag = (nSelCorrvTag / nEvtsv);
    accCorrv = (nSelCorrv / nEvtsv);

    accErrCorrvFSR = (accCorrvFSR * sqrt((nSelCorrVarv) / (nSelCorrvFSR * nSelCorrvFSR) + 1. / nEvtsv));
    accErrCorrvMC = (accCorrvMC * sqrt((nSelCorrVarv) / (nSelCorrvMC * nSelCorrvMC) + 1. / nEvtsv));
    accErrCorrvBkg = (accCorrvBkg * sqrt((nSelCorrVarv) / (nSelCorrvBkg * nSelCorrvBkg) + 1. / nEvtsv));
    accErrCorrvTag = (accCorrvTag * sqrt((nSelCorrVarv) / (nSelCorrvTag * nSelCorrvTag) + 1. / nEvtsv));
    //accErrCorrv.push_back(accCorrv[ifile] * sqrt((nSelCorrVarv[ifile]) / (nSelCorrv[ifile] * nSelCorrv[ifile]) + 1. / nEvtsv[ifile]));

    //--------------------------------------------------------------------------------------------------------------
    // Output
    //==============================================================================================================
    // Print full set for efficiency calculations
    char masterOutput[600];
    // just start printing....
    //for (uint ifile = 0; ifile < fnamev.size(); ++ifile) { // go through info per file
    for (uint ifile = 0; ifile < 1; ++ifile) {
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
    cout << "*" << endl;
    cout << "* SUMMARY" << endl;
    cout << "*--------------------------------------------------" << endl;
    if (charge == 0)
        cout << " W -> e nu" << endl;
    if (charge == -1)
        cout << " W- -> e nu" << endl;
    if (charge == 1)
        cout << " W+ -> e nu" << endl;
    cout << "  pT > " << PT_CUT << endl;
    cout << "  |eta| < " << ETA_CUT << endl;
    cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
    cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
    cout << endl;

    double acc_sys = 0.;
    acc_sys += TMath::Power(accErrCorrv_pos / accCorrv, 2.0);
    acc_sys += TMath::Power(accErrCorrv_neg / accCorrv, 2.0);
    acc_sys += TMath::Power(accCorrvFSR / accCorrv - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvMC / accCorrv - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvBkg / accCorrv - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvTag / accCorrv - 1.0, 2.0);
    acc_sys = TMath::Sqrt(acc_sys);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (uint ifile = 0; ifile < 1; ++ifile) {
        cout << "    *** Acceptance ***" << endl;
        cout << "     barrel: " << setw(12) << nSelBv << " / " << nEvtsv << " = " << accBv << " +/- " << accErrBv;
        cout << "  ==eff corr==> " << accBCorrv << " +/- " << accErrBCorrv << endl;
        cout << "     endcap: " << setw(12) << nSelEv << " / " << nEvtsv << " = " << accEv << " +/- " << accErrEv;
        cout << "  ==eff corr==> " << accECorrv << " +/- " << accErrECorrv << endl;
        cout << "      total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv;
        cout << "  ==eff corr==> " << accCorrv << " +/- " << accErrCorrv << endl;
        cout << "  ==total efficiency==> " << setw(4) << accCorrv / accv << endl;
        cout << "          stat: " << 100 * accErrCorrv / accCorrv << endl;
        cout << "          stat pos: " << 100 * accErrCorrv_pos / accCorrv << endl;
        cout << "          stat neg: " << 100 * accErrCorrv_neg / accCorrv << endl;
        cout << "          FSR unc: " << 100 * (accCorrvFSR / accCorrv - 1.0) << endl;
        cout << "           MC unc: " << 100 * (accCorrvMC / accCorrv - 1.0)  << endl;
        cout << "          Bkg unc: " << 100 * (accCorrvBkg / accCorrv - 1.0) << endl;
        cout << "          Tag unc: " << 100 * (accCorrvTag / accCorrv - 1.0) << endl;
        cout << "          Total: " << 100 * acc_sys << endl;
        cout << endl;

        cout << endl << endl;
        cout << " Prefire Correction " << nSelPfire / nSelCorrv << " + " << fabs(nSelPfireUp - nSelPfire) / nSelPfire << " - " << fabs(nSelPfireDown - nSelPfire) / nSelPfire << endl;
        cout << " Prefire ECAL Correction " << nSelPfireEcal / nSelCorrv << " + " << fabs(nSelPfireEcalUp - nSelPfireEcal) / nSelPfireEcal << " - " << fabs(nSelPfireEcalDown - nSelPfireEcal) / nSelPfireEcal << endl;
        cout << " Prefire Muon Correction " << nSelPfireMuon / nSelCorrv << " + " << fabs(nSelPfireMuonUp - nSelPfireMuon) / nSelPfireMuon << " - " << fabs(nSelPfireMuonDown - nSelPfireMuon) / nSelPfireMuon << endl;
        cout << " Prefire Photon Correction " << nSelPfirePhoton / nSelCorrv << " + " << fabs(nSelPfirePhotonUp - nSelPfirePhoton) / nSelPfirePhoton << " - " << fabs(nSelPfirePhotonDown - nSelPfirePhoton) / nSelPfirePhoton << endl;
        cout << " Prefire Jet Correction " << nSelPfireJet / nSelCorrv << " + " << fabs(nSelPfireJetUp - nSelPfireJet) / nSelPfireJet << " - " << fabs(nSelPfireJetDown - nSelPfireJet) / nSelPfireJet << endl;
    }

    char txtfname[500];
    sprintf(txtfname, "%s/sel.txt", outputDir.Data());
    ofstream txtfile;
    txtfile.open(txtfname);
    txtfile << "*" << endl;
    txtfile << "* SUMMARY" << endl;
    txtfile << "*--------------------------------------------------" << endl;
    if (charge == 0)
        txtfile << " W -> e nu" << endl;
    if (charge == -1)
        txtfile << " W- -> e nu" << endl;
    if (charge == 1)
        txtfile << " W+ -> e nu" << endl;
    txtfile << "  pT > " << PT_CUT << endl;
    txtfile << "  |eta| < " << ETA_CUT << endl;
    txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
    txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
    txtfile << endl;

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (uint ifile = 0; ifile < 1; ++ifile) {
        txtfile << "    *** Acceptance ***" << endl;
        txtfile << "     barrel: " << setw(12) << nSelBv << " / " << nEvtsv << " = " << accBv << " +/- " << accErrBv;
        txtfile << "  ==eff corr==> " << accBCorrv << " +/- " << accErrBCorrv << endl;
        txtfile << "     endcap: " << setw(12) << nSelEv << " / " << nEvtsv << " = " << accEv << " +/- " << accErrEv;
        txtfile << "  ==eff corr==> " << accECorrv << " +/- " << accErrECorrv << endl;
        txtfile << "      total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv;
        txtfile << "  ==eff corr==> " << accCorrv << " +/- " << accErrCorrv << endl;
        txtfile << "  ==total efficiency==> " << setw(4) << accCorrv / accv << endl;
        txtfile << "     SF corrected: " << accCorrv << " +/- " << accErrCorrv << endl;
        txtfile << "     SF corrected Pos: " << accCorrv << " +/- " << accErrCorrv_pos << endl;
        txtfile << "     SF corrected Neg: " << accCorrv << " +/- " << accErrCorrv_neg << endl;
        txtfile << "          FSR unc: " << accCorrvFSR << " +/- " << accErrCorrvFSR << endl;
        txtfile << "  ==pct diff (FSR) ==> " << 100 * (accCorrvFSR - accCorrv) / accCorrv << endl;
        txtfile << "           MC unc: " << accCorrvMC << " +/- " << accErrCorrvMC << endl;
        txtfile << "  ==pct diff (MC) ==> " << 100 * (accCorrvMC - accCorrv) / accCorrv << endl;
        txtfile << "          Bkg unc: " << accCorrvBkg << " +/- " << accErrCorrvBkg << endl;
        txtfile << "  ==pct diff (bkg) ==> " << 100 * (accCorrvBkg - accCorrv) / accCorrv << endl;
        txtfile << endl;
    }
    txtfile.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sel_nums_only.txt", outputDir.Data());
    ofstream txtfile2;
    txtfile2.open(txtfname);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (uint ifile = 0; ifile < 1; ++ifile) {
        txtfile2 << "uncorrected: " << accv << endl;
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

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (uint ifile = 0; ifile < 1; ++ifile) {
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

    gBenchmark->Show("computeAccSelWe");
}
