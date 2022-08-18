//================================================================================================
//
// Compute W->munu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties need to be checked//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TLorentzVector.h"
#include <TBenchmark.h> // class to track macro running statistics
#include <TClonesArray.h> // ROOT array class
#include <TFile.h> // file handle class
#include <TH1D.h> // histogram class
#include <TROOT.h> // access to gROOT, entry point to ROOT system
#include <TSystem.h> // interface to OS
#include <TTree.h> // class to access ntuples
#include <fstream> // functions for file I/O
#include <iomanip> // functions to format standard I/O
#include <iostream> // standard I/O
#include <sstream> // class for parsing strings
#include <string> // C++ string class
#include <vector> // STL vector class

#include <numeric>

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

#include "MitEwk13TeV/Utils/CSample.hh" // helper class to handle samples
#include "MitEwk13TeV/Utils/ConfParse.hh" // input conf file parser
#include "MitEwk13TeV/Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "MitEwk13TeV/Utils/MyTools.hh" // various helper functions

// helper class to handle efficiency tables
#include "MitEwk13TeV/Utils/AppEffSF.cc"
#include "MitEwk13TeV/Utils/CEffUser2D.hh"
#endif

//=== MAIN MACRO =================================================================================================

void computeAccSelWm(const TString conf, // input file
    const TString inputDir,
    const TString outputDir, // output directory
    const TString outputName,
    const Int_t charge, // 0 = inclusive, +1 = W+, -1 = W-
    const Int_t doPU,
    const TString sysFileSIT, // condense these into 1 file per type of eff (pos & neg into 1 file)
    const TString sysFileSta,
    const bool is13TeV = 1)
{
    gBenchmark->Start("computeAccSelWm");

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================
    const Double_t mu_MASS = 0.1057;

    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;

    const Double_t ETA_BARREL = 1.442;
    const Double_t ETA_ENDCAP = 1.442;
    const Double_t VETO_PT = 10;
    const Double_t VETO_ETA = 2.4;

    const Int_t BOSON_ID = 24;
    const Int_t LEPTON_ID = 13;

    const int NBptSta = 1;
    const float ptrangeSta[NBptSta + 1] = {25., 100000.};

    const int NBetaSta = 18;
    const float etarangeSta[NBetaSta + 1] = {-2.4, -2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,-0.15,0,0.15,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};

    const int NBptSIT = 4;
    const float ptrangeSIT[NBptSIT + 1] = { 25., 30, 35, 40, 100000.};

    const int NBeta = 12;
    const float etarange[NBeta + 1] = { -2.4, -2.1, -1.6, -1.2, -0.9, -0.3, 0, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4 };

    const int NBptHLT = 12;
    const float ptrangeHLT[NBptHLT + 1] = { 25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 100000.0 };

    AppEffSF effs(inputDir);
    effs.loadHLT("MuHLTEff_aMCxPythia", "Positive", "Negative");
    effs.loadSel("MuSITEff_aMCxPythia", "Combined", "Combined");
    effs.loadSta("MuStaEff_aMCxPythia", "Combined", "Combined");
    //effs.loadHLT("MuHLTEff_aMCxPythia", "Positive", "Negative");
    //effs.loadSel("MuSISEff_aMCxPythia", "Combined", "Combined");
    //effs.loadSta("MuTrkEff_aMCxPythia", "Combined", "Combined");
    effs.loadUncSel(sysFileSIT);
    effs.loadUncSta(sysFileSta);

    TH2D* hSelErr_pos = new TH2D("hSelErr_pos", "", NBeta, etarange, NBptSIT, ptrangeSIT);
    TH2D* hSelErr_neg = new TH2D("hSelErr_neg", "", NBeta, etarange, NBptSIT, ptrangeSIT);

    TH2D* hStaErr_pos = new TH2D("hStaErr_pos", "", NBetaSta, etarangeSta, NBptSta, ptrangeSta);
    TH2D* hStaErr_neg = new TH2D("hStaErr_neg", "", NBetaSta, etarangeSta, NBptSta, ptrangeSta);

    TH2D* hHLTErr_pos = new TH2D("hHLTErr_pos", "", NBeta, etarange, NBptHLT, ptrangeHLT);
    TH2D* hHLTErr_neg = new TH2D("hHLTErr_neg", "", NBeta, etarange, NBptHLT, ptrangeHLT);

    //// load pileup reweighting file
    //TFile* f_rw = TFile::Open("../Tools/pileup_rw_76X.root", "read");
    //TH1D* h_rw = (TH1D*)f_rw->Get("h_rw_golden");

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
    // Data structures to store info from TTrees
    baconhep::TEventInfo* info = new baconhep::TEventInfo();
    baconhep::TGenEventInfo* gen = new baconhep::TGenEventInfo();
    TClonesArray* genPartArr = new TClonesArray("baconhep::TGenParticle");
    TClonesArray* muonArr = new TClonesArray("baconhep::TMuon");
    TClonesArray* vertexArr = new TClonesArray("baconhep::TVertex");

    TFile* infile = 0;
    TTree* eventTree = 0;

    // Variables to store acceptances and uncertainties
    Double_t nEvtsv = 0, nSelv = 0, nSelBv = 0, nSelEv = 0, nAllv = 0;
    Double_t accv = 0, accBv = 0, accEv = 0;
    Double_t accErrv = 0, accErrBv = 0, accErrEv = 0;
    Double_t nSelCorrv = 0, nSelBCorrv = 0, nSelECorrv = 0;
    Double_t nSelRaw = 0, nSelCorrToHLT = 0, nSelCorrToSel = 0, nSelCorrToTrk = 0;
    Double_t nSelCorrVarv = 0, nSelBCorrVarv = 0, nSelECorrVarv = 0, nSelCorrVarvPos = 0, nSelCorrVarvNeg = 0;
    Double_t accCorrv = 0, accBCorrv = 0, accECorrv = 0;
    Double_t accErrCorrv = 0, accErrBCorrv = 0, accErrECorrv = 0, accErrCorrv_pos = 0, accErrCorrv_neg = 0;

    Double_t nSelCorrvFSR = 0, nSelCorrvMC = 0, nSelCorrvBkg = 0, nSelCorrvTag = 0; //, nSelCorrvStat;
    Double_t nSelCorrvFSR_I = 0, nSelCorrvMC_I = 0, nSelCorrvBkg_I = 0, nSelCorrvTag_I = 0; //, nSelCorrvStat_I;
    Double_t nSelCorrvFSR_S = 0, nSelCorrvMC_S = 0, nSelCorrvBkg_S = 0, nSelCorrvTag_S = 0; //, nSelCorrvStat_S;
    Double_t accCorrvFSR = 0, accCorrvMC = 0, accCorrvBkg = 0, accCorrvTag = 0; //, accCorrvStat;
    Double_t accCorrvFSR_I = 0, accCorrvMC_I = 0, accCorrvBkg_I = 0, accCorrvTag_I = 0; //, accCorrvStat_I;
    Double_t accCorrvFSR_S = 0, accCorrvMC_S = 0, accCorrvBkg_S = 0, accCorrvTag_S = 0; //, accCorrvStat_S;
    Double_t accErrCorrvFSR = 0, accErrCorrvMC = 0, accErrCorrvBkg = 0, accErrCorrvTag = 0; //, accErrCorrvStat;

    const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

    CSample* samp = 0;

    //
    // loop over the processes
    //
    for (UInt_t isamp = 0; isamp < samplev.size(); ++isamp) {

        //
        // loop through files
        //
        CSample* samp = samplev[isamp];
        const UInt_t nfiles = samp->fnamev.size();

        for (UInt_t ifile = 0; ifile < nfiles; ifile++) {

            // Read input file and get the TTrees
            cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... ";
            cout.flush();
            infile = TFile::Open(samp->fnamev[ifile]);

            eventTree = (TTree*)infile->Get("Events");
            assert(eventTree);
            eventTree->SetBranchAddress("Info", &info);
            TBranch* infoBr = eventTree->GetBranch("Info");
            eventTree->SetBranchAddress("GenEvtInfo", &gen);
            TBranch* genBr = eventTree->GetBranch("GenEvtInfo");
            eventTree->SetBranchAddress("GenParticle", &genPartArr);
            TBranch* genPartBr = eventTree->GetBranch("GenParticle");
            eventTree->SetBranchAddress("Muon", &muonArr);
            TBranch* muonBr = eventTree->GetBranch("Muon");
            eventTree->SetBranchAddress("PV", &vertexArr);
            TBranch* vertexBr = eventTree->GetBranch("PV");

            //
            // loop over events
            //
            double frac = 0.05; // fraction of events to be used for calculation
            if (isamp == 0) frac = 0.05;
            if (isamp == 1) frac = 0.10;
            if (isamp == 2) frac = 0.30;
            if (is13TeV) frac = 0.30;
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
                genBr->GetEntry(ientry);
                infoBr->GetEntry(ientry);
                genPartArr->Clear();
                genPartBr->GetEntry(ientry);

                if (fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                    continue;
                if (charge == -1 && toolbox::flavor(genPartArr, BOSON_ID) != LEPTON_ID)
                    continue;
                if (charge == 1 && toolbox::flavor(genPartArr, BOSON_ID) != -LEPTON_ID)
                    continue;
                if (charge == 0 && fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                    continue;

                Double_t weight = (gen->weight > 0 ? 1 : -1) * samp->xsecv[ifile] / nWgtSum;
                nAllv += weight;

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

                double genMuon_pt = 0;
                double genMuon_eta = 0;
                double genMuon_phi = 0;
                if( ((toolbox::flavor(genPartArr, BOSON_ID) == LEPTON_ID) && (glepq1 < 0)) || ((toolbox::flavor(genPartArr, BOSON_ID) == -LEPTON_ID) && (glepq1 > 0)) ) {
                    genMuon_pt = glep1->Pt();
                    genMuon_eta = glep1->Eta();
                    genMuon_phi = glep1->Phi();
                } else {
                    genMuon_pt = glep2->Pt();
                    genMuon_eta = glep2->Eta();
                    genMuon_phi = glep2->Phi();
                }

                //if (genMuon_pt < PT_CUT)
                //    continue;
                //if (fabs(genMuon_eta) > ETA_CUT)
                //    continue;

                // TLorentzVector tvec=*glep1+*glep2;
                // TLorentzVector* genV=new TLorentzVector(0,0,0,0);
                // genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
                // genVPt   = tvec.Pt();
                // genVPhi  = tvec.Phi();
                // genVy    = tvec.Rapidity();
                // double genVMass = tvec.M();
                //double mtgen = sqrt(2.0 * (glep1->Pt()) * (glep2->Pt()) * (1.0 - cos(toolbox::deltaPhi(glep1->Phi(), glep2->Phi()))));
                // if(mtgen > 40) continue;
                //if (mtgen < 40)
                //    continue;
                // cout << "mass " << genVMass <<  "   mt " << mt << endl;
                // cout << "-- out 2 " << endl;
                vertexArr->Clear();
                vertexBr->GetEntry(ientry);
                double npv = vertexArr->GetEntries();
                //if (doPU > 0)
                //    weight *= h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

                nEvtsv += weight;

                TLorentzVector vMu(0, 0, 0, 0);
                // trigger requirement
                // if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
                //if (!isMuonTrigger(triggerMenu, info->triggerBits, kFALSE, is13TeV))
                //    continue;

                // good vertex requirement
                if (!(info->hasGoodPV))
                    continue;
                // cout << "-- out 3 " << endl;
                
                muonArr->Clear();
                muonBr->GetEntry(ientry);
                Int_t nLooseLep = 0;
                const baconhep::TMuon* goodMuon = 0;
                Bool_t passSel = kFALSE;

                for (Int_t i = 0; i < muonArr->GetEntriesFast(); i++) {
                    const baconhep::TMuon* mu = (baconhep::TMuon*)((*muonArr)[i]);

                    if (fabs(mu->eta) > VETO_ETA)
                        continue; // loose lepton |eta| cut
                    if (mu->pt < VETO_PT)
                        continue; // loose lepton pT cut
                    if (passMuonLooseID(mu))
                        nLooseLep++; // loose lepton selection
                    if (nLooseLep > 1) { // extra lepton veto
                        passSel = kFALSE;
                        break;
                    }

                    if (fabs(mu->eta) > ETA_CUT)
                        continue; // lepton |eta| cut
                    if (mu->pt < PT_CUT)
                        continue; // lepton pT cut

                    //cout << "mu pt " << mu->pt << " eta " << mu->eta << " phi " << mu->phi << " gen mu pt " << genMuon_pt << " eta" << genMuon_eta << " phi " << genMuon_phi << endl;
                    //if( toolbox::deltaR(mu->eta, mu->phi, genMuon_eta, genMuon_phi) > 0.5)
                    //    continue;

                    if (!passMuonID(mu))
                        continue; // lepton selection

                    if (!isMuonTriggerObj(triggerMenu, mu->hltMatchBits, kFALSE, is13TeV))
                        continue;

                    if (!isMuonTrigger(triggerMenu, info->triggerBits, kFALSE, is13TeV))
                        continue;


                    if (charge != 0 && mu->q != charge)
                        continue; // check charge (if necessary)
                    // cout << "-- out 5 " << endl;

                    double mtreco = sqrt(2.0 * (mu->pt) * (info->pfMETC) * (1.0 - cos(toolbox::deltaPhi(mu->phi, info->pfMETCphi))));

                    if (mtreco < 20)
                        continue;
                    // if(mtreco > 40 && mtreco < 140) continue;
                    // cout << "-- out 6 " << endl;
                    passSel = kTRUE;
                    goodMuon = mu;
                    vMu.SetPtEtaPhiM(mu->pt, mu->eta, mu->phi, mu_MASS);
                }

                if (passSel) {

                    /******** We have a W candidate! HURRAY! ********/

                    Bool_t isBarrel = (fabs(goodMuon->eta) < ETA_BARREL) ? kTRUE : kFALSE;

                    // data/MC scale factor corrections
                    Double_t effdata, effmc, emTag;
                    Double_t edFSR, edMC, edBkg, edTag; //, edStat;
                    Double_t corr = 1;
                    Double_t corrFSR = 1, corrMC = 1, corrBkg = 1, corrTag = 1; //, corrStat=1;
                    Double_t corrFSR_I = 1, corrMC_I = 1, corrBkg_I = 1, corrTag_I = 1; //, corrStat_I=1;
                    Double_t corrFSR_S = 1, corrMC_S = 1, corrBkg_S = 1, corrTag_S = 1; //, corrStat_S=1;

                    effdata = 1;
                    effmc = 1;
                    emTag = 1;
                    edFSR = 1;
                    edMC = 1;
                    edBkg = 1;
                    edTag = 1; //edStat=1;

                    int q = goodMuon->q;

                    corr = effs.fullCorrections(&vMu, q);
                    // corr = effs.dataOnly(&vMu,q);
                    vector<double> uncs_sta = effs.getUncSta(&vMu, q);
                    vector<double> uncs_sit = effs.getUncSel(&vMu, q);

                    corrFSR *= uncs_sta[0] * uncs_sit[0] * effs.computeHLTSF(&vMu, q); // alternate fsr model
                    corrMC *= uncs_sta[1] * uncs_sit[1] * effs.computeHLTSF(&vMu, q); // alternate mc gen model
                    corrBkg *= uncs_sta[2] * uncs_sit[2] * effs.computeHLTSF(&vMu, q); // alternate bkg model
                    corrTag *= uncs_sta[3] * uncs_sit[3] * effs.computeHLTSF(&vMu, q); // alternate bkg model
                    // corr *= effdata/effmc; // orig

                    corrFSR_I *= uncs_sit[0] * effs.computeHLTSF(&vMu, q) * effs.computeStaSF(&vMu, q); // alternate fsr model
                    corrMC_I *= uncs_sit[1] * effs.computeHLTSF(&vMu, q) * effs.computeStaSF(&vMu, q); // alternate mc gen model
                    corrBkg_I *= uncs_sit[2] * effs.computeHLTSF(&vMu, q) * effs.computeStaSF(&vMu, q); // alternate bkg model
                    corrTag_I *= uncs_sit[3] * effs.computeHLTSF(&vMu, q) * effs.computeStaSF(&vMu, q); // alternate bkg model

                    corrFSR_S *= uncs_sta[0] * effs.computeHLTSF(&vMu, q) * effs.computeSelSF(&vMu, q); // alternate fsr model
                    corrMC_S *= uncs_sta[1] * effs.computeHLTSF(&vMu, q) * effs.computeSelSF(&vMu, q); // alternate mc gen model
                    corrBkg_S *= uncs_sta[2] * effs.computeHLTSF(&vMu, q) * effs.computeSelSF(&vMu, q); // alternate bkg model
                    corrTag_S *= uncs_sta[3] * effs.computeHLTSF(&vMu, q) * effs.computeSelSF(&vMu, q); // alternate bkg model

                    double var = 0.;
                    // var += effs.statUncSta(&l1, q) + effs.statUncSta(&l2, q2);
                    var += effs.statUncSta(&vMu, q, hStaErr_pos, hStaErr_neg, fabs(weight) * corr, true);
                    var += effs.statUncSel(&vMu, q, hSelErr_pos, hSelErr_neg, fabs(weight) * corr, true);
                    var += effs.statUncHLT(&vMu, q, hHLTErr_pos, hHLTErr_neg, fabs(weight) * corr);

                    nSelv += weight;
                    nSelCorrvFSR += weight * corrFSR;
                    nSelCorrvFSR_I += weight * corrFSR_I;
                    nSelCorrvFSR_S += weight * corrFSR_S;
                    nSelCorrvMC += weight * corrMC;
                    nSelCorrvMC_I += weight * corrMC_I;
                    nSelCorrvMC_S += weight * corrMC_S;
                    nSelCorrvBkg += weight * corrBkg;
                    nSelCorrvBkg_I += weight * corrBkg_I;
                    nSelCorrvBkg_S += weight * corrBkg_S;
                    nSelCorrvTag += weight * corrTag;
                    nSelCorrvTag_I += weight * corrTag_I;
                    nSelCorrvTag_S += weight * corrTag_S;
                    // /*nSelCorrvStat[ifile]+=weight*corrStat; */nSelCorrvStat_I[ifile]+=weight*corrStat_I; nSelCorrvStat_S[ifile]+=weight*corrStat_S;

                    nSelCorrv += weight * corr;
                    // nSelCorrVarv[ifile]+=weight*weight*corr*corr;
                    
                    double mu_pt = vMu.Pt();
                    double mu_eta = vMu.Eta();

                    double eff_mc_hlt = 1.0;
                    double eff_mc_sta = 1.0;
                    double eff_mc_sit = 1.0;
                    if (q > 0) {
                        eff_mc_hlt = effs.hlt.mcPos.getEff(mu_eta, mu_pt);
                        eff_mc_sta = effs.sta.mcPos.getEff(mu_eta, mu_pt);
                        eff_mc_sit = effs.sel.mcPos.getEff(mu_eta, mu_pt);
                        effs.effSel(&vMu, q, 0);
                    } else {
                        eff_mc_hlt = effs.hlt.mcNeg.getEff(mu_eta, mu_pt);
                        eff_mc_sta = effs.sta.mcNeg.getEff(mu_eta, mu_pt);
                        eff_mc_sit = effs.sel.mcNeg.getEff(mu_eta, mu_pt);
                    }
                    eff_mc_hlt = 1.0;
                    nSelCorrToHLT += weight / eff_mc_hlt;
                    nSelCorrToSel += weight / (eff_mc_hlt * eff_mc_sit);
                    nSelCorrToTrk += weight / (eff_mc_hlt * eff_mc_sit * eff_mc_sta);
                    nSelRaw += weight;

                    if (isBarrel) {
                        nSelBv += weight;
                        nSelBCorrv += weight * corr;
                        nSelBCorrVarv += weight * weight * corr * corr;
                    } else {
                        nSelEv += weight;
                        nSelECorrv += weight * corr;
                        nSelECorrVarv += weight * weight * corr * corr;
                    }
                }
            }

            //delete infile;
            //infile = 0, eventTree = 0;
            std::cout << " nSelCorrVarv " << nSelCorrVarv << " nSelCorrVarvPos " << nSelCorrVarvPos << " nSelCorrVarvNeg " << nSelCorrVarvNeg << std::endl;
        }
        //delete info;
        //delete gen;
        //delete muonArr;
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
            // std::cout << "hlt pos " << var << std::endl;
        }
    }

    for (Int_t iy = 0; iy <= hSelErr_pos->GetNbinsY() + 1; iy++) {
        for (Int_t ix = 0; ix <= hSelErr_pos->GetNbinsX() + 1; ix++) {
            Double_t err = hSelErr_pos->GetBinContent(ix, iy);
            var += err * err;
            var_pos += err * err;
            err = hSelErr_neg->GetBinContent(ix, iy);
            var += err * err;
            var_neg += err * err;
            // std::cout << "sel pos " << var << std::endl;
        }
    }

    for (Int_t iy = 0; iy <= hStaErr_pos->GetNbinsY() + 1; iy++) {
        for (Int_t ix = 0; ix <= hStaErr_pos->GetNbinsX() + 1; ix++) {
            Double_t err = hStaErr_pos->GetBinContent(ix, iy);
            var += err * err;
            var_pos += err * err;
            err = hStaErr_neg->GetBinContent(ix, iy);
            var += err * err;
            var_neg += err * err;
        }
    }
    nSelCorrVarv += var;
    nSelCorrVarvPos += var_pos;
    nSelCorrVarvNeg += var_neg;

    // compute acceptances
    accv = (nSelv / nEvtsv);
    accErrv = (sqrt(accv * (1. + accv) / nEvtsv));
    accBv = (nSelBv / nEvtsv);
    accErrBv = (sqrt(accBv * (1. + accBv) / nEvtsv));
    accEv = (nSelEv / nEvtsv);
    accErrEv = (sqrt(accEv * (1. + accEv) / nEvtsv));

    accCorrv = (nSelCorrv / nEvtsv);
    std::cout << "variation  " << nSelCorrVarv / (nSelCorrv * nSelCorrv) << std::endl;
    //accErrCorrv = (accCorrv * sqrt(nSelCorrVarv / (nSelCorrv * nSelCorrv) + 1. / nEvtsv));
    //accErrCorrv_pos = (accCorrv * sqrt(nSelCorrVarvPos / (nSelCorrv * nSelCorrv) + 1. / nEvtsv));
    //accErrCorrv_neg = (accCorrv * sqrt(nSelCorrVarvNeg / (nSelCorrv * nSelCorrv) + 1. / nEvtsv));
    accErrCorrv = (accCorrv * sqrt(nSelCorrVarv) / nSelCorrv );
    accErrCorrv_pos = (accCorrv * sqrt(nSelCorrVarvPos) / nSelCorrv );
    accErrCorrv_neg = (accCorrv * sqrt(nSelCorrVarvNeg) / nSelCorrv );

    accBCorrv = (nSelBCorrv / nEvtsv);
    accErrBCorrv = (accBCorrv * sqrt(nSelBCorrVarv / nSelBCorrv / nSelBCorrv + 1. / nEvtsv));
    accECorrv = (nSelECorrv / nEvtsv);
    accErrECorrv = (accECorrv * sqrt(nSelECorrVarv / (nSelECorrv * nSelECorrv) + 1. / nEvtsv));

    accCorrvFSR = (nSelCorrvFSR / nEvtsv);
    accCorrvMC = (nSelCorrvMC / nEvtsv);
    accCorrvBkg = (nSelCorrvBkg / nEvtsv);
    accCorrvTag = (nSelCorrvTag / nEvtsv);

    accCorrvFSR_I = (nSelCorrvFSR_I / nEvtsv);
    accCorrvMC_I = (nSelCorrvMC_I / nEvtsv);
    accCorrvBkg_I = (nSelCorrvBkg_I / nEvtsv);
    accCorrvTag_I = (nSelCorrvTag_I / nEvtsv);

    accCorrvFSR_S = (nSelCorrvFSR_S / nEvtsv);
    accCorrvMC_S = (nSelCorrvMC_S / nEvtsv);
    accCorrvBkg_S = (nSelCorrvBkg_S / nEvtsv);
    accCorrvTag_S = (nSelCorrvTag_S / nEvtsv);

    accErrCorrvFSR = (accCorrvFSR * sqrt((nSelCorrVarv) / (nSelCorrvFSR * nSelCorrvFSR) + 1. / nEvtsv));
    accErrCorrvMC = (accCorrvMC * sqrt((nSelCorrVarv) / (nSelCorrvMC * nSelCorrvMC) + 1. / nEvtsv));
    accErrCorrvBkg = (accCorrvBkg * sqrt((nSelCorrVarv) / (nSelCorrvBkg * nSelCorrvBkg) + 1. / nEvtsv));
    accErrCorrvTag = (accCorrvTag * sqrt((nSelCorrVarv) / (nSelCorrvTag * nSelCorrvTag) + 1. / nEvtsv));

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
                << " " << accCorrvFSR_I << endl;
        txtfile << "sta_fsr"
                << " " << accCorrvFSR_S << endl;
        txtfile << "sel_mc"
                << " " << accCorrvMC_I << endl;
        txtfile << "sta_mc"
                << " " << accCorrvMC_S << endl;
        txtfile << "sel_bkg"
                << " " << accCorrvBkg_I << endl;
        txtfile << "sta_bkg"
                << " " << accCorrvBkg_S << endl;
        txtfile << "sel_tagpt"
                << " " << accCorrvTag_I << endl;
        txtfile << "sta_tagpt"
                << " " << accCorrvTag_S << endl;
        txtfile.close();
    }

    cout << "*" << endl;
    cout << "* SUMMARY" << endl;
    cout << "*--------------------------------------------------" << endl;
    if (charge == 0)
        cout << " W -> mu nu" << endl;
    if (charge == -1)
        cout << " W- -> mu nu" << endl;
    if (charge == 1)
        cout << " W+ -> mu nu" << endl;
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
        cout << "   ================================================" << endl;
        cout << "     *** counting *** " << endl;
        cout << "     After Gen selection: " << nAllv << endl;
        cout << "     after Gen lepton selection " << nEvtsv << endl;
        cout << "     **** " << endl;
        cout << endl;

        // 2D binned scale factor
        cout << "     **** " << endl;
        cout << "     Raw Counting: " << nSelRaw << endl;
        cout << "     SF corrected: " << nSelCorrv << endl;
        cout << "     HLT corrected: " << nSelCorrToHLT << endl;
        cout << "     Sit corrected: " << nSelCorrToSel << endl;
        cout << "     Trk corrected: " << nSelCorrToTrk << endl;
        cout << endl;

        cout << "    *** Acceptance ***" << endl;
        cout << "     barrel: " << setw(12) << nSelBv << " / " << nEvtsv << " = " << accBv << " +/- " << accErrBv;
        cout << "  ==eff corr==> " << accBCorrv << " +/- " << accErrBCorrv << endl;
        cout << "     endcap: " << setw(12) << nSelEv << " / " << nEvtsv << " = " << accEv << " +/- " << accErrEv;
        cout << "  ==eff corr==> " << accECorrv << " +/- " << accErrECorrv << endl;
        cout << "      total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv;
        cout << "  ==eff corr==> " << accCorrv << " +/- " << accErrCorrv << endl;
        cout << "  ==total efficiency==> " << setw(4) << accCorrv / accv << endl;
        cout << "     SF corrected: " << accCorrv << " +/- " << accErrCorrv << endl;
        cout << "          stat: " << 100 * accErrCorrv / accCorrv << endl;
        cout << "          stat pos: " << accErrCorrv_pos / accCorrv << endl;
        cout << "          stat neg: " << accErrCorrv_neg / accCorrv << endl;
        cout << "          FSR unc: " << accCorrvFSR / accCorrv - 1.0 << " / Sel: " << accCorrvFSR_I << " / Sta: " << accCorrvFSR_S << endl;
        cout << "           MC unc: " << accCorrvMC / accCorrv - 1.0 << " / Sel: " << accCorrvMC_I << " / Sta: " << accCorrvMC_S << endl;
        cout << "          Bkg unc: " << accCorrvBkg / accCorrv - 1.0 << " / Sel: " << accCorrvBkg_I << " / Sta: " << accCorrvBkg_S << endl;
        cout << "          Tag unc: " << accCorrvTag / accCorrv - 1.0 << " / Sel: " << accCorrvTag_I << " / Sta: " << accCorrvTag_S << endl;
        cout << "          Total: " << acc_sys << endl;
        cout << "Acc (FSR/MC/Bkg/Tag): " << accCorrvFSR << ", " << accCorrvMC << ", " << accCorrvBkg << ", " << accCorrvTag << endl;
        cout << "  fraction passing gen cut: " << nEvtsv << " / " << nAllv << " = " << nEvtsv / nAllv << endl;
        cout << endl;
    }

    char txtfname[500];
    sprintf(txtfname, "%s/sel.txt", outputDir.Data());
    ofstream txtfile;
    txtfile.open(txtfname);
    txtfile << "*" << endl;
    txtfile << "* SUMMARY" << endl;
    txtfile << "*--------------------------------------------------" << endl;
    if (charge == 0)
        txtfile << " W -> mu nu" << endl;
    if (charge == -1)
        txtfile << " W- -> mu nu" << endl;
    if (charge == 1)
        txtfile << " W+ -> mu nu" << endl;
    txtfile << "  pT > " << PT_CUT << endl;
    txtfile << "  |eta| < " << ETA_CUT << endl;
    txtfile << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
    txtfile << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
    txtfile << endl;

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (uint ifile = 0; ifile < 1; ++ifile) {
        txtfile << "   ================================================" << endl;
        txtfile << "    *** Acceptance ***" << endl;
        txtfile << "     barrel: " << setw(12) << nSelBv << " / " << nEvtsv << " = " << accBv << " +/- " << accErrBv;
        txtfile << "  ==eff corr==> " << accBCorrv << " +/- " << accErrBCorrv << endl;
        txtfile << "     endcap: " << setw(12) << nSelEv << " / " << nEvtsv << " = " << accEv << " +/- " << accErrEv;
        txtfile << "  ==eff corr==> " << accECorrv << " +/- " << accErrECorrv << endl;
        txtfile << "      total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv;
        txtfile << "  ==eff corr==> " << accCorrv << " +/- " << accErrCorrv << endl;
        txtfile << "  ==total efficiency==> " << setw(4) << accCorrv / accv << endl;
        txtfile << "     SF corrected: " << accCorrv << " +/- " << accErrCorrv << endl;
        txtfile << "          FSR unc: " << accCorrvFSR << " / Sel: " << accCorrvFSR_I << " / Sta: " << accCorrvFSR_S << endl;
        txtfile << "           MC unc: " << accCorrvMC << " / Sel: " << accCorrvMC_I << " / Sta: " << accCorrvMC_S << endl;
        txtfile << "          Bkg unc: " << accCorrvBkg << " / Sel: " << accCorrvBkg_I << " / Sta: " << accCorrvBkg_S << endl;
        txtfile << "          Tag unc: " << accCorrvTag << " / Sel: " << accCorrvTag_I << " / Sta: " << accCorrvTag_S << endl;
        txtfile << "Acc (FSR/MC/Bkg/Tag): " << accCorrvFSR << ", " << accCorrvMC << ", " << accCorrvBkg << ", " << accCorrvTag << endl;
        txtfile << "  fraction passing gen cut: " << nEvtsv << " / " << nAllv << " = " << nEvtsv / nAllv << endl;

        txtfile << endl;
    }
    txtfile.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sel_nums_only.txt", outputDir.Data());
    ofstream txtfile2;
    txtfile2.open(txtfname);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) { //accv[ifile]
    for (uint ifile = 0; ifile < 1; ++ifile) {
        txtfile2 << "uncorrected: " << accv << endl;
        txtfile2 << accCorrv << " " << accErrCorrv << endl;
        txtfile2 << accCorrvFSR << " " << accCorrvFSR_I << " " << accCorrvFSR_S << endl;
        txtfile2 << accCorrvMC << " " << accCorrvMC_I << " " << accCorrvMC_S << endl;
        txtfile2 << accCorrvBkg << " " << accCorrvBkg_I << " " << accCorrvBkg_S << endl;
        txtfile2 << accCorrvTag << " " << accCorrvTag_I << " " << accCorrvTag_S << endl;
        // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

        txtfile2 << endl;
    }
    txtfile2.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sit_unc.txt", outputDir.Data());
    ofstream txtfile3;
    txtfile3.open(txtfname);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (uint ifile = 0; ifile < 1; ++ifile) {
        txtfile3 << accCorrv << endl;
        txtfile3 << accCorrvFSR_I << endl;
        txtfile3 << accCorrvMC_I << endl;
        txtfile3 << accCorrvBkg_I << endl;
        txtfile3 << accCorrvTag_I << endl;
        // txtfile << accCorrvFSR  << ", " << accCorrvMC << ", " << accCorrvBkg << ", " << accCorrvTag << endl;

        // txtfile3 << endl;
    }
    txtfile3.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sta_unc.txt", outputDir.Data());
    ofstream txtfile4;
    txtfile4.open(txtfname);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (uint ifile = 0; ifile < 1; ++ifile) {
        txtfile4 << accCorrv << endl;
        txtfile4 << accCorrvFSR_S << endl;
        txtfile4 << accCorrvMC_S << endl;
        txtfile4 << accCorrvBkg_S << endl;
        txtfile4 << accCorrvTag_S << endl;
        // txtfile << accCorrvFSR  << ", " << accCorrvMC << ", " << accCorrvBkg << ", " << accCorrvTag << endl;

        // txtfile4 << endl;
    }
    txtfile4.close();

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;

    gBenchmark->Show("computeAccSelWm");
}
