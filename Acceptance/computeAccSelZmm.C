//================================================================================================
//
// Compute Z->mumu acceptance at full selection level
//
//  * outputs results summary text file
//
//  [!!!] propagation of efficiency scale factor uncertainties no yet implemented
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
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
// #include <math.h>                  // mathematics
#include "TLorentzVector.h" // 4-vector class

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
#include "MitEwk13TeV/RochesterCorr/RoccoR.cc"

// helper class to handle efficiency tables
#include "MitEwk13TeV/Utils/AppEffSF.cc"
#include "MitEwk13TeV/Utils/CEffUser1D.hh"
#include "MitEwk13TeV/Utils/CEffUser2D.hh"
#endif

//=== MAIN MACRO =================================================================================================

void computeAccSelZmm(const TString conf, // input file
    const TString inputDir,
    const TString outputDir, // output directory
    const TString outputName,
    const Int_t doPU,
    const TString sysFileSIT, // condense these into 1 file per type of eff (pos & neg into 1 file)
    const TString sysFileSta,
    const bool is13TeV = 1)
{
    gBenchmark->Start("computeAccSelZmm");

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================

    const Double_t MASS_LOW = 60;
    const Double_t MASS_HIGH = 120;
    const Double_t PT_CUT = 25;
    const Double_t PT_CUT_PRESELECT = 20;
    const Double_t ETA_CUT = 2.4;
    const Double_t MUON_MASS = 0.105658369;

    const Int_t BOSON_ID = 23;
    const Int_t LEPTON_ID = 13;

    const int NBptSta = 1;
    // max pt should be a large number, 
    // so that the overflow bin can be correlated with the max bin
    const float ptrangeSta[NBptSta + 1] = {25., 100000.};

    const int NBetaSta = 18;
    const float etarangeSta[NBetaSta + 1] = {-2.4, -2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,-0.15,0,0.15,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};

    const int NBptSIT = 4;
    const float ptrangeSIT[NBptSIT + 1] = { 25., 30, 35, 40, 100000.};

    const int NBeta = 12;
    const float etarange[NBeta + 1] = { -2.4, -2.1, -1.6, -1.2, -0.9, -0.3, 0, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4 };

    const int NBptHLT = 12;
    const float ptrangeHLT[NBptHLT + 1] = { 25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 1000000 };

    AppEffSF effs(inputDir);
    effs.loadHLT("MuHLTEff_aMCxPythia", "Positive", "Negative");
    effs.loadSel("MuSITEff_aMCxPythia", "Combined", "Combined");
    effs.loadSta("MuStaEff_aMCxPythia", "Combined", "Combined");
    effs.loadUncSel(sysFileSIT);
    effs.loadUncSta(sysFileSta);

    // load pileup reweighting file
    //TFile* f_rw = TFile::Open("../Tools/puWeights_76x.root", "read");
    //TH1D* h_rw = (TH1D*)f_rw->Get("puWeights");

    RoccoR rc("/uscms/home/yfeng/nobackup/WpT/CMSSW_9_4_19/src/MitEwk13TeV/RochesterCorr/RoccoR2017.txt");

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

    TH2D* hSelErr_pos = new TH2D("hSelErr_pos", "", NBeta, etarange, NBptSIT, ptrangeSIT);
    TH2D* hSelErr_neg = new TH2D("hSelErr_neg", "", NBeta, etarange, NBptSIT, ptrangeSIT);

    TH2D* hStaErr_pos = new TH2D("hStaErr_pos", "", NBeta, etarange, NBptSta, ptrangeSta);
    TH2D* hStaErr_neg = new TH2D("hStaErr_neg", "", NBeta, etarange, NBptSta, ptrangeSta);

    TH2D* hHLTErr_pos = new TH2D("hHLTErr_pos", "", NBeta, etarange, NBptHLT, ptrangeHLT);
    TH2D* hHLTErr_neg = new TH2D("hHLTErr_neg", "", NBeta, etarange, NBptHLT, ptrangeHLT);

    // Data structures to store info from TTrees
    baconhep::TEventInfo* info = new baconhep::TEventInfo();
    baconhep::TGenEventInfo* gen = new baconhep::TGenEventInfo();
    TClonesArray* genPartArr = new TClonesArray("baconhep::TGenParticle");
    TClonesArray* muonArr = new TClonesArray("baconhep::TMuon");
    TClonesArray* vertexArr = new TClonesArray("baconhep::TVertex");

    TFile* infile = 0;
    TTree* eventTree = 0;

    // Variables to store acceptances and uncertainties (per input file)
    vector<Double_t> nEvtsv, nSelv, nAllv;
    vector<Double_t> nEvtsv_inFiducial, nSelv_inFiducial, nSelCorrv_inFiducial, nSelCorrv_MC, nSelCorrv_lepup, nSelCorrv_lepdn, nSelCorrv_leps;
    vector<Double_t> nEvtsv_pv, nEvtsv_reco, nEvtsv_sta, nEvtsv_sel, nEvtsv_trig;
    vector<Double_t> nSelCorrv, nSelCorrVarv, nSelCorrVarv_pos, nSelCorrVarv_neg;
    vector<Double_t> accv, accCorrv;
    vector<Double_t> accErrv, accErrCorrv, accErrCorrv_pos, accErrCorrv_neg;
    vector<Double_t> nSelCorrvFSR, nSelCorrvMC, nSelCorrvBkg, nSelCorrvTag; //, nSelCorrvStat;
    vector<Double_t> nSelCorrvFSR_I, nSelCorrvMC_I, nSelCorrvBkg_I, nSelCorrvTag_I; //, nSelCorrvStat_I;
    vector<Double_t> nSelCorrvFSR_S, nSelCorrvMC_S, nSelCorrvBkg_S, nSelCorrvTag_S; //, nSelCorrvStat_S;
    vector<Double_t> nSelCorrVarvFSR, nSelCorrVarvMC, nSelCorrVarvBkg, nSelCorrVarvTag; //, nSelCorrVarvStat;
    vector<Double_t> nSelCorrv_MC_trg, nSelCorrv_MC_sel, nSelCorrv_MC_sit; // for closure test
    vector<Double_t> accCorrvFSR, accCorrvMC, accCorrvBkg, accCorrvTag; //, accCorrvStat;
    vector<Double_t> accCorrvFSR_I, accCorrvMC_I, accCorrvBkg_I, accCorrvTag_I; //, accCorrvStat_I;
    vector<Double_t> accCorrvFSR_S, accCorrvMC_S, accCorrvBkg_S, accCorrvTag_S; //, accCorrvStat_S;
    vector<Double_t> accErrCorrvFSR, accErrCorrvMC, accErrCorrvBkg, accErrCorrvTag; //, accErrCorrvStat;
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
        eventTree->SetBranchAddress("Muon", &muonArr);
        TBranch* muonBr = eventTree->GetBranch("Muon");
        eventTree->SetBranchAddress("PV", &vertexArr);
        TBranch* vertexBr = eventTree->GetBranch("PV");

        nAllv.push_back(0);
        nEvtsv.push_back(0);
        nEvtsv_inFiducial.push_back(0);
        nEvtsv_pv.push_back(0);
        nEvtsv_reco.push_back(0);
        nEvtsv_sta.push_back(0);
        nEvtsv_sel.push_back(0);
        nEvtsv_trig.push_back(0);
        nSelv.push_back(0);
        nSelv_inFiducial.push_back(0);
        nSelCorrv_inFiducial.push_back(0);
        nSelCorrv_MC.push_back(0);
        nSelCorrv_lepup.push_back(0);
        nSelCorrv_lepdn.push_back(0);
        nSelCorrv_leps.push_back(0);
        nSelCorrv.push_back(0);
        nSelCorrVarv.push_back(0);
        nSelCorrVarv_pos.push_back(0);
        nSelCorrVarv_neg.push_back(0);
        nSelCorrvFSR.push_back(0);
        nSelCorrVarvFSR.push_back(0);
        nSelCorrvMC.push_back(0);
        nSelCorrVarvMC.push_back(0);
        nSelCorrvBkg.push_back(0);
        nSelCorrVarvBkg.push_back(0);
        nSelCorrvTag.push_back(0);
        nSelCorrVarvTag.push_back(0);
        nSelCorrvFSR_I.push_back(0);
        nSelCorrvFSR_S.push_back(0);
        nSelCorrvMC_I.push_back(0);
        nSelCorrvMC_S.push_back(0);
        nSelCorrvBkg_I.push_back(0);
        nSelCorrvBkg_S.push_back(0);
        nSelCorrvTag_I.push_back(0);
        nSelCorrvTag_S.push_back(0);

        nSelCorrv_MC_trg.push_back(0);
        nSelCorrv_MC_sel.push_back(0);
        nSelCorrv_MC_sit.push_back(0);
        //
        // loop over events
        //
        double frac = 0.10;
        if (is13TeV) 
            frac = 0.30;
        for (UInt_t ientry = 0; ientry < (uint)(frac * eventTree->GetEntries()); ientry++) {
            if (ientry % 100000 == 0)
                cout << "Processing event " << ientry << ". " << (double)ientry / (double)eventTree->GetEntries() * 100 << " percent done with this file." << endl;
            genBr->GetEntry(ientry);
            genPartArr->Clear();
            genPartBr->GetEntry(ientry);
            infoBr->GetEntry(ientry);

            Int_t glepq1 = -99;
            Int_t glepq2 = -99;

            if (fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                continue;

            //nAllv[ifile] += gen->weight;
            nAllv[ifile] += (gen->weight > 0) ? 1: -1;

            TLorentzVector* vec = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector* lep1 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector* lep2 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector* lep3 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector* lep4 = new TLorentzVector(0, 0, 0, 0);
            toolbox::fillGenBorn(genPartArr, BOSON_ID, vec, lep1, lep2, lep3, lep4);

            bool doDressed = true;

            if (doDressed) {
                //cout << "do dressed lepton " << endl;
                TLorentzVector* gph = new TLorentzVector(0, 0, 0, 0);
                for (Int_t i = 0; i < genPartArr->GetEntriesFast(); i++) {
                    const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*)((*genPartArr)[i]);
                    if (fabs(genloop->pdgId) != 22)
                        continue;
                    gph->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                    if (toolbox::deltaR(gph->Eta(), gph->Phi(), lep3->Eta(), lep3->Phi()) < 0.1)
                        lep3->operator+=(*gph);
                    if (toolbox::deltaR(gph->Eta(), gph->Phi(), lep4->Eta(), lep4->Phi()) < 0.1)
                        lep4->operator+=(*gph);
                }
                delete gph;
            }

            double genmass =vec->M();
            //if (ientry%1==0)
            //    cout << "mass " << genmass << endl;

            // if((vec->M()<MASS_LOW || vec->M()>MASS_HIGH)) continue;
            vertexArr->Clear();
            vertexBr->GetEntry(ientry);
            double npv = vertexArr->GetEntries();
            //Double_t weight = gen->weight;
            Double_t weight = (gen->weight > 0) ? 1: -1;
            //if (doPU > 0)
            //    weight *= h_rw->GetBinContent(h_rw->FindBin(info->nPUmean));

            nEvtsv[ifile] += weight;

            // trigger requirement
            if (!isMuonTrigger(triggerMenu, info->triggerBits, kFALSE, is13TeV))
                continue;

            // good vertex requirement
            if (!(info->hasGoodPV))
                continue;
            nEvtsv_pv[ifile] += weight;

            muonArr->Clear();
            muonBr->GetEntry(ientry);

            for (Int_t i1 = 0; i1 < muonArr->GetEntriesFast(); i1++) {
                double weight_l1 = weight;
                double weight_l1_lepup = weight;
                double weight_l1_lepdn = weight;
                double weight_l1_leps = weight;

                const baconhep::TMuon* mu1 = (baconhep::TMuon*)((*muonArr)[i1]);

                if (mu1->pt < PT_CUT_PRESELECT)
                    continue; // lepton pT cut
                if (fabs(mu1->eta) > ETA_CUT)
                    continue; // lepton |eta| cut
                if (!passMuonID(mu1))
                    continue; // lepton selection

                int q1 = mu1->q;

                TLorentzVector vMu1(0, 0, 0, 0);
                vMu1.SetPtEtaPhiM(mu1->pt, mu1->eta, mu1->phi, MUON_MASS);
                TLorentzVector vMu1u, vMu1d;
                vMu1u = vMu1;
                vMu1d = vMu1;

                // rochester correction
                float genMuonPt1 = -1;
                if (lep3 && toolbox::deltaR(mu1->eta, mu1->phi, lep3->Eta(), lep3->Phi()) < 0.3) {
                    genMuonPt1 = lep3->Pt();
                } else if (lep4 && toolbox::deltaR(mu1->eta, mu1->phi, lep4->Eta(), lep4->Phi()) < 0.3) {
                    genMuonPt1 = lep4->Pt();
                }
                double mcSF1 = 1;
                double deltaMcSF1 = 0.;
                if (genMuonPt1 > 0) {
                    mcSF1 = rc.kSpreadMC(q1, vMu1.Pt(), vMu1.Eta(), vMu1.Phi(), genMuonPt1);
                    deltaMcSF1 = rc.kSpreadMCerror(q1, vMu1.Pt(), vMu1.Eta(), vMu1.Phi(), genMuonPt1);
                } else {
                    double rand = gRandom->Uniform(1);
                    mcSF1 = rc.kSmearMC(q1, vMu1.Pt(), vMu1.Eta(), vMu1.Phi(), mu1->nTkLayers, rand);
                    deltaMcSF1 = rc.kSmearMCerror(q1, vMu1.Pt(), vMu1.Eta(), vMu1.Phi(), mu1->nTkLayers, rand);
                }
                mcSF1 = 1.0;
                vMu1 *= mcSF1;
                vMu1u *= mcSF1 * (1 + deltaMcSF1);
                vMu1d *= mcSF1 * (1 - deltaMcSF1);
                if (vMu1.Pt() < PT_CUT)
                    continue; // lepton pT cut
                // do not use continue option to evaluate the lepton 
                // energy scale effects
                double weight_lepup = weight;
                double weight_lepdn = weight;

                weight_l1_lepup *= (vMu1u.Pt() >= PT_CUT);
                weight_l1_lepdn *= (vMu1d.Pt() >= PT_CUT);

                double rsmear = gRandom->Gaus(0,1);
                // check the smearing effect on the Z->mumu counts
                TLorentzVector vMu1s;
                vMu1s.SetPtEtaPhiM( vMu1.Pt() + rsmear * 0.3, vMu1.Eta(), vMu1.Phi(), vMu1.M() );
                weight_l1_leps *= (vMu1s.Pt() >= PT_CUT);

                for (Int_t i2 = i1 + 1; (i2 < muonArr->GetEntriesFast()); i2++) {

                    double weight_l2 = weight_l1;
                    double weight_l2_lepup = weight_l1_lepup;
                    double weight_l2_lepdn = weight_l1_lepdn;
                    double weight_l2_leps = weight_l1_leps;

                    const baconhep::TMuon* mu2 = (baconhep::TMuon*)((*muonArr)[i2]);

                    if (mu1->q == mu2->q)
                        continue; // opposite charge requirement
                    if (mu2->pt < PT_CUT_PRESELECT)
                        continue; // lepton pT cut
                    if (fabs(mu2->eta) > ETA_CUT)
                        continue; // lepton |eta| cut
                    if (!passMuonID(mu2))
                        continue; // lepton selection

                    int q2 = mu2->q;

                    TLorentzVector vMu2(0, 0, 0, 0);
                    vMu2.SetPtEtaPhiM(mu2->pt, mu2->eta, mu2->phi, MUON_MASS);
                    TLorentzVector vMu2u, vMu2d;
                    vMu2u = vMu2;
                    vMu2d = vMu2;

                    // rochester correction for mu2
                    float genMuonPt2 = -1;
                    if (lep3 && toolbox::deltaR(mu2->eta, mu2->phi, lep3->Eta(), lep3->Phi()) < 0.3) {
                        genMuonPt2 = lep3->Pt();
                    } else if (lep4 && toolbox::deltaR(mu2->eta, mu2->phi, lep4->Eta(), lep4->Phi()) < 0.3) {
                        genMuonPt2 = lep4->Pt();
                    }
                    double mcSF2 = 1;
                    double deltaMcSF2 = 0.;
                    if (genMuonPt2 > 0) {
                        mcSF2 = rc.kSpreadMC(q2, vMu2.Pt(), vMu2.Eta(), vMu2.Phi(), genMuonPt2);
                        deltaMcSF2 = rc.kSpreadMCerror(q2, vMu2.Pt(), vMu2.Eta(), vMu2.Phi(), genMuonPt2);
                    } else {
                        double rand = gRandom->Uniform(1);
                        mcSF2 = rc.kSmearMC(q2, vMu2.Pt(), vMu2.Eta(), vMu2.Phi(), mu2->nTkLayers, rand);
                        deltaMcSF2 = rc.kSmearMCerror(q2, vMu2.Pt(), vMu2.Eta(), vMu2.Phi(), mu2->nTkLayers, rand);
                    }

                    //std::cout << "mcSF1 " << mcSF1 << " deltaMcSF1 " << deltaMcSF1 << " mcSF2 " << mcSF2 << " deltaMcSF2 " << deltaMcSF2 << std::endl;
                    mcSF2 = 1.0;
                    vMu2 *= mcSF2;
                    vMu2u *= mcSF2 * (1 + deltaMcSF2);
                    vMu2d *= mcSF2 * (1 - deltaMcSF2);
                    if (vMu2.Pt() < PT_CUT)
                        continue; // lepton pT cut
                    // do not use continue option to evaluate the lepton
                    // energy scale effects
                    weight_l2_lepup *= (vMu2u.Pt() >= PT_CUT);
                    weight_l2_lepdn *= (vMu2d.Pt() >= PT_CUT);

                    double rsmear = gRandom->Gaus(0,1);
                    TLorentzVector vMu2s;
                    vMu2s.SetPtEtaPhiM( vMu2.Pt() + rsmear * 0.3, vMu2.Eta(), vMu2.Phi(), vMu2.M() );
                    weight_l2_leps *= (vMu2s.Pt() >= PT_CUT);


                    // mass window
                    TLorentzVector vDilep = vMu1 + vMu2;
                    if ((vDilep.M() < MASS_LOW) || (vDilep.M() > MASS_HIGH))
                        continue;

                    TLorentzVector vDilepu = vMu1u + vMu2u;
                    weight_l2_lepup *= ((vDilepu.M() >= MASS_LOW) && (vDilepu.M() <= MASS_HIGH));
                    TLorentzVector vDilepd = vMu1d + vMu2d;
                    weight_l2_lepdn *= ((vDilepd.M() >= MASS_LOW) && (vDilepd.M() <= MASS_HIGH));

                    TLorentzVector vDileps = vMu1s + vMu2s;
                    weight_l2_leps *= ((vDileps.M() >= MASS_LOW) && (vDileps.M() <= MASS_HIGH));

                    //std::cout << "mu1 pt " << vMu1.Pt() << " up " << vMu1u.Pt() << " down " << vMu1d.Pt() << std::endl;
                    //std::cout << "mu2 pt " << vMu2.Pt() << " up " << vMu2u.Pt() << " down " << vMu2d.Pt() << std::endl;
                    //std::cout << "invariant mass " << vDilep.M() << " up " << vDilepu.M() << " down " << vDilepd.M() << std::endl;

                    // gen reco matching
                    bool isMatched = false;
                    Bool_t match1 = false;
                    Bool_t match2 = false;
                    if (lep3 && lep4 && lep3->Pt() > 0 && lep4->Pt() > 0) {
                        match1 = ((toolbox::deltaR(vMu1.Eta(), vMu1.Phi(), lep3->Eta(), lep3->Phi()) < 0.3) || (toolbox::deltaR(vMu1.Eta(), vMu1.Phi(), lep4->Eta(), lep4->Phi()) < 0.3));
                        match2 = ((toolbox::deltaR(vMu2.Eta(), vMu2.Phi(), lep3->Eta(), lep3->Phi()) < 0.3) || (toolbox::deltaR(vMu2.Eta(), vMu2.Phi(), lep4->Eta(), lep4->Phi()) < 0.3));
                    }
                    if (match1 && match2 && q1!=q2 )
                        isMatched = true;

                    // trigger match
                    if (!isMuonTriggerObj(triggerMenu, mu1->hltMatchBits, kFALSE, is13TeV) && !isMuonTriggerObj(triggerMenu, mu2->hltMatchBits, kFALSE, is13TeV))
                        continue;

                    /******** We have a Z candidate! HURRAY! ********/
                    Double_t corrFSR = 1, corrMC = 1, corrBkg = 1, corrTag = 1; //, corrStat=1;
                    Double_t corrFSR_I = 1, corrMC_I = 1, corrBkg_I = 1, corrTag_I = 1; //, corrStat_I=1;
                    Double_t corrFSR_S = 1, corrMC_S = 1, corrBkg_S = 1, corrTag_S = 1; //, corrStat_S=1;


                    double corr = effs.fullCorrections(&vMu1, q1, &vMu2, q2);
                    // corr = effs.dataOnly(&vMu1,q1,&vMu2,q2);
                    vector<double> uncs_sta = effs.getUncSta(&vMu1, q1, &vMu2, q2);
                    vector<double> uncs_sit = effs.getUncSel(&vMu1, q1, &vMu2, q2);

                    corrFSR *= uncs_sta[0] * uncs_sit[0] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2); // alternate fsr model
                    corrMC *= uncs_sta[1] * uncs_sit[1] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2); // alternate mc gen model
                    corrBkg *= uncs_sta[2] * uncs_sit[2] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2); // alternate bkg model
                    corrTag *= uncs_sta[3] * uncs_sit[3] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2); // alternate tag-pt cut

                    corrFSR_I *= uncs_sit[0] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeStaSF(&vMu1, q1, &vMu2, q2);
                    corrMC_I *= uncs_sit[1] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeStaSF(&vMu1, q1, &vMu2, q2);
                    corrBkg_I *= uncs_sit[2] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeStaSF(&vMu1, q1, &vMu2, q2);
                    corrTag_I *= uncs_sit[3] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeStaSF(&vMu1, q1, &vMu2, q2);

                    corrFSR_S *= uncs_sta[0] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeSelSF(&vMu1, q1, &vMu2, q2);
                    corrMC_S *= uncs_sta[1] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeSelSF(&vMu1, q1, &vMu2, q2);
                    corrBkg_S *= uncs_sta[2] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeSelSF(&vMu1, q1, &vMu2, q2);
                    corrTag_S *= uncs_sta[3] * effs.computeHLTSF(&vMu1, q1, &vMu2, q2) * effs.computeSelSF(&vMu1, q1, &vMu2, q2);

                    double var = 0.;
                    // var += effs.statUncSta(&l1, q1) + effs.statUncSta(&l2, q2);
                    // should it be fabs(weight) or weight?
                    // i.e., should the negative weight sample contribute positively to the stat unc, or negatively?
                    var += effs.statUncSta(&vMu1, q1, hStaErr_pos, hStaErr_neg, fabs(weight_l2) * corr, true);
                    var += effs.statUncSta(&vMu2, q2, hStaErr_pos, hStaErr_neg, fabs(weight_l2) * corr, true);
                    var += effs.statUncSel(&vMu1, q1, hSelErr_pos, hSelErr_neg, fabs(weight_l2) * corr, true);
                    var += effs.statUncSel(&vMu2, q2, hSelErr_pos, hSelErr_neg, fabs(weight_l2) * corr, true);
                    var += effs.statUncHLT(&vMu1, q1, hHLTErr_pos, hHLTErr_neg, fabs(weight_l2) * corr);
                    var += effs.statUncHLT(&vMu2, q2, hHLTErr_pos, hHLTErr_neg, fabs(weight_l2) * corr);

                    // std::cout << info->evtNum << " " << corr << " " << std::endl;
                    nSelv[ifile] += weight_l2;
                    nSelCorrvFSR[ifile] += weight_l2 * corrFSR;
                    nSelCorrvFSR_I[ifile] += weight_l2 * corrFSR_I;
                    nSelCorrvFSR_S[ifile] += weight_l2 * corrFSR_S;
                    nSelCorrvMC[ifile] += weight_l2 * corrMC;
                    nSelCorrvMC_I[ifile] += weight_l2 * corrMC_I;
                    nSelCorrvMC_S[ifile] += weight_l2 * corrMC_S;
                    nSelCorrvBkg[ifile] += weight_l2 * corrBkg;
                    nSelCorrvBkg_I[ifile] += weight_l2 * corrBkg_I;
                    nSelCorrvBkg_S[ifile] += weight_l2 * corrBkg_S;
                    nSelCorrvTag[ifile] += weight_l2 * corrTag;
                    nSelCorrvTag_I[ifile] += weight_l2 * corrTag_I;
                    nSelCorrvTag_S[ifile] += weight_l2 * corrTag_S;

                    // get the MC efficiencies
                    double effMC = effs.efficiencyOnly(&vMu1, q1, &vMu2, q2, 0);
                    nSelCorrv_MC[ifile] += weight_l2 / effMC;
                    nSelCorrv[ifile] += weight_l2 * corr;
                    nSelCorrVarvFSR[ifile] += weight_l2 * weight_l2 * corrFSR * corrFSR;
                    nSelCorrVarvMC[ifile] += weight_l2 * weight_l2 * corrMC * corrMC;
                    nSelCorrVarvBkg[ifile] += weight_l2 * weight_l2 * corrBkg * corrBkg;
                    nSelCorrVarvTag[ifile] += weight_l2 * weight_l2 * corrTag * corrTag;
                    // nSelCorrVarv[ifile]+=weight*weight*corr*corr;
                    //
                    nSelCorrv_lepup[ifile] += weight_l2_lepup * corr;
                    nSelCorrv_lepdn[ifile] += weight_l2_lepdn * corr;
                    nSelCorrv_leps[ifile]  += weight_l2_leps * corr;
                    
                    // for closure test
                    double effMC_sta = effs.effSta(&vMu1, q1, &vMu2, q2, 0);
                    double effMC_sel = effs.effSel(&vMu1, q1, &vMu2, q2, 0);
                    double effMC_trg = effs.effHLT(&vMu1, q1, &vMu2, q2, 0);
                    nSelCorrv_MC_trg[ifile] += weight_l2 / effMC_trg * isMatched;
                    nSelCorrv_MC_sel[ifile] += weight_l2 / (effMC_trg * effMC_sel) * isMatched;
                    nSelCorrv_MC_sit[ifile] += weight_l2 / (effMC_trg * effMC_sel * effMC_sta) * isMatched;
                }
            }

            delete vec;
            delete lep1;
            delete lep2;
            delete lep3;
            delete lep4;
        }

        // std::cout << "nSelCorrVarv[ifile]  " <<  nSelCorrVarv[ifile] << std::endl;
        // std::cout << "var" << std::endl;
        Double_t var = 0;
        Double_t var_pos = 0;
        Double_t var_neg = 0;

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

        std::cout << "var " << var << " pos " << var_pos << " neg " << var_neg << std::endl;

        nSelCorrVarv[ifile] += var;
        nSelCorrVarv_pos[ifile] += var_pos;
        nSelCorrVarv_neg[ifile] += var_neg;
        // compute acceptances
        accv.push_back(nSelv[ifile] / nEvtsv[ifile]);
        accErrv.push_back(accv[ifile] * sqrt((1. + accv[ifile]) / nEvtsv[ifile]));

        accCorrvFSR.push_back(nSelCorrvFSR[ifile] / nEvtsv[ifile]);
        accCorrvMC.push_back(nSelCorrvMC[ifile] / nEvtsv[ifile]);
        accCorrvBkg.push_back(nSelCorrvBkg[ifile] / nEvtsv[ifile]);
        accCorrvTag.push_back(nSelCorrvTag[ifile] / nEvtsv[ifile]);

        accCorrvFSR_I.push_back(nSelCorrvFSR_I[ifile] / nEvtsv[ifile]);
        accCorrvMC_I.push_back(nSelCorrvMC_I[ifile] / nEvtsv[ifile]);
        accCorrvBkg_I.push_back(nSelCorrvBkg_I[ifile] / nEvtsv[ifile]);
        accCorrvTag_I.push_back(nSelCorrvTag_I[ifile] / nEvtsv[ifile]);

        accCorrvFSR_S.push_back(nSelCorrvFSR_S[ifile] / nEvtsv[ifile]);
        accCorrvMC_S.push_back(nSelCorrvMC_S[ifile] / nEvtsv[ifile]);
        accCorrvBkg_S.push_back(nSelCorrvBkg_S[ifile] / nEvtsv[ifile]);
        accCorrvTag_S.push_back(nSelCorrvTag_S[ifile] / nEvtsv[ifile]);

        std::cout << "statistical variation " << sqrt((nSelCorrVarv[ifile]) / (nSelCorrv[ifile] * nSelCorrv[ifile])) << std::endl;
        accErrCorrvFSR.push_back(accCorrvFSR[ifile] * sqrt((nSelCorrVarv[ifile]) / (nSelCorrvFSR[ifile] * nSelCorrvFSR[ifile]) + 1. / nEvtsv[ifile]));
        accErrCorrvMC.push_back(accCorrvMC[ifile] * sqrt((nSelCorrVarv[ifile]) / (nSelCorrvMC[ifile] * nSelCorrvMC[ifile]) + 1. / nEvtsv[ifile]));
        accErrCorrvBkg.push_back(accCorrvBkg[ifile] * sqrt((nSelCorrVarv[ifile]) / (nSelCorrvBkg[ifile] * nSelCorrvBkg[ifile]) + 1. / nEvtsv[ifile]));
        accErrCorrvTag.push_back(accCorrvTag[ifile] * sqrt((nSelCorrVarv[ifile]) / (nSelCorrvTag[ifile] * nSelCorrvTag[ifile]) + 1. / nEvtsv[ifile]));

        accCorrv.push_back(nSelCorrv[ifile] / nEvtsv[ifile]);
        //accErrCorrv.push_back(accCorrv[ifile] * sqrt((nSelCorrVarv[ifile]) / (nSelCorrv[ifile] * nSelCorrv[ifile]) + 1. / nEvtsv[ifile]));
        accErrCorrv.push_back(accCorrv[ifile] * sqrt(nSelCorrVarv[ifile]) / nSelCorrv[ifile] );

        //accErrCorrv_pos.push_back(accCorrv[ifile] * sqrt((nSelCorrVarv_pos[ifile]) / (nSelCorrv[ifile] * nSelCorrv[ifile]) + 1. / nEvtsv[ifile]));
        accErrCorrv_pos.push_back(accCorrv[ifile] * sqrt(nSelCorrVarv_pos[ifile]) / nSelCorrv[ifile] );
        //accErrCorrv_neg.push_back(accCorrv[ifile] * sqrt((nSelCorrVarv_neg[ifile]) / (nSelCorrv[ifile] * nSelCorrv[ifile]) + 1. / nEvtsv[ifile]));
        accErrCorrv_neg.push_back( accCorrv[ifile] * sqrt(nSelCorrVarv_neg[ifile]) / nSelCorrv[ifile] );

        delete infile;
        infile = 0, eventTree = 0;
    }

    delete info;
    delete gen;
    delete muonArr;

    // Print full set for efficiency calculations
    char masterOutput[600];
    // just start printing....
    //for (uint ifile = 0; ifile < fnamev.size(); ++ifile) { // go through info per file
    for (uint ifile = 0; ifile < 1; ++ifile) { // go through info per file
        sprintf(masterOutput, "%s/%s.txt", outputDir.Data(), outputName.Data());
        ofstream txtfile;
        txtfile.open(masterOutput);
        txtfile << "acc " << accCorrv[ifile] << endl;
        txtfile << "acc_stat " << accCorrv[ifile] + accErrCorrv[ifile] << endl;
        txtfile << "sel_fsr"
                << " " << accCorrvFSR_I[ifile] << endl;
        txtfile << "sta_fsr"
                << " " << accCorrvFSR_S[ifile] << endl;
        txtfile << "sel_mc"
                << " " << accCorrvMC_I[ifile] << endl;
        txtfile << "sta_mc"
                << " " << accCorrvMC_S[ifile] << endl;
        txtfile << "sel_bkg"
                << " " << accCorrvBkg_I[ifile] << endl;
        txtfile << "sta_bkg"
                << " " << accCorrvBkg_S[ifile] << endl;
        txtfile << "sel_tagpt"
                << " " << accCorrvTag_I[ifile] << endl;
        txtfile << "sta_tagpt"
                << " " << accCorrvTag_S[ifile] << endl;
        txtfile.close();
    }

    std::cout << "print output" << std::endl;
    //--------------------------------------------------------------------------------------------------------------
    // Output
    //==============================================================================================================
    cout << "*" << endl;
    cout << "* SUMMARY" << endl;
    cout << "*--------------------------------------------------" << endl;
    cout << " Z -> mu mu" << endl;
    cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
    cout << "  pT > " << PT_CUT << endl;
    cout << "  |eta| < " << ETA_CUT << endl;
    cout << endl;

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    double acc_sys = 0.;
    int ifile = 0;
    acc_sys += TMath::Power(accErrCorrv_pos[ifile] / accCorrv[ifile], 2.0);
    acc_sys += TMath::Power(accErrCorrv_neg[ifile] / accCorrv[ifile], 2.0);
    acc_sys += TMath::Power(accCorrvFSR[ifile] / accCorrv[ifile] - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvMC[ifile] / accCorrv[ifile] - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvBkg[ifile] / accCorrv[ifile] - 1.0, 2.0);
    acc_sys += TMath::Power(accCorrvTag[ifile] / accCorrv[ifile] - 1.0, 2.0);
    acc_sys = TMath::Sqrt(acc_sys);

    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        cout << "     After PV selection: " << nEvtsv_pv[ifile] << endl;
        cout << "     After RecoMuon selection: " << nEvtsv_reco[ifile] << endl;
        cout << "     After standalone muon requirement " << nEvtsv_sta[ifile] << endl;
        cout << "     After Muon selection: " << nEvtsv_sel[ifile] << endl;
        cout << "     After trig requirement: " << nEvtsv_trig[ifile] << endl;
        cout << endl;
        cout << "     efficiency correction with MC trig " << nSelCorrv_MC_trg[ifile] << endl;
        cout << "     efficiency correction with MC trig * sel " << nSelCorrv_MC_sel[ifile] << endl;
        cout << "     efficiency correction with MC trig * sel * sit " << nSelCorrv_MC_sit[ifile] << endl;
        cout << endl;
        cout << "    *** Acceptance ***" << endl;
        cout << "          nominal: " << setw(12) << nSelv[ifile] << " / " << nEvtsv[ifile] << " = " << accv[ifile] << " +/- " << accErrv[ifile]  << ", after correction " << nSelCorrv_MC[ifile] << endl;
        cout << "          In the fiducial region " << setw(12) << nSelv_inFiducial[ifile] << " / " << nEvtsv_inFiducial[ifile] << ", after correction " << nSelCorrv_inFiducial[ifile] << endl;
        cout << "          with lepton momentum scale variation up " << setw(12) << nSelCorrv_lepup[ifile] << " and down " << nSelCorrv_lepdn[ifile] << " central " << nSelCorrv[ifile] << " smearing " << nSelCorrv_leps[ifile] << endl;
        cout << "     SF corrected: " << accCorrv[ifile] << " +/- " << accErrCorrv[ifile] << endl;
        cout << "     SF corrected pos: " << accCorrv[ifile] << " +/- " << accErrCorrv_pos[ifile] << endl;
        cout << "     SF corrected neg: " << accCorrv[ifile] << " +/- " << accErrCorrv_neg[ifile] << endl;
        cout << "  ==total efficiency==> " << setw(4) << accCorrv[ifile] / accv[ifile] << endl;
        cout << "          stat: " << 100 * accErrCorrv[ifile] / accCorrv[ifile] << endl;
        cout << "          stat pos: " << accErrCorrv_pos[ifile] / accCorrv[ifile] << endl;
        cout << "          stat neg: " << accErrCorrv_neg[ifile] / accCorrv[ifile] << endl;
        
        cout << "          FSR unc: " << accCorrvFSR[ifile] / accCorrv[ifile] - 1.0 << " / Sel: " << accCorrvFSR_I[ifile] << " / Sta: " << accCorrvFSR_S[ifile] << endl;
        cout << "           MC unc: " << accCorrvMC[ifile] / accCorrv[ifile] - 1.0 << " / Sel: " << accCorrvMC_I[ifile] << " / Sta: " << accCorrvMC_S[ifile] << endl;
        cout << "          Bkg unc: " << accCorrvBkg[ifile] / accCorrv[ifile] - 1.0 << " / Sel: " << accCorrvBkg_I[ifile] << " / Sta: " << accCorrvBkg_S[ifile] << endl;
        cout << "          Tag unc: " << accCorrvTag[ifile] / accCorrv[ifile] - 1.0 << " / Sel: " << accCorrvTag_I[ifile] << " / Sta: " << accCorrvTag_S[ifile] << endl;
        // cout << "         Stat unc: " << accCorrvStat[ifile] << " / Sel: " << accCorrvStat_I[ifile]<< " / Sta: " << accCorrvStat_S[ifile] << endl;
        cout << "          Total: " << acc_sys << endl;
        cout << "  fraction passing gen cut: " << nEvtsv[ifile] << " / " << nAllv[ifile] << " = " << nEvtsv[ifile] / nAllv[ifile] << endl;
        cout << endl;
    }

    char txtfname[500];
    sprintf(txtfname, "%s/binned.txt", outputDir.Data());
    ofstream txtfile;
    txtfile.open(txtfname);
    txtfile << "*" << endl;
    txtfile << "* SUMMARY" << endl;
    txtfile << "*--------------------------------------------------" << endl;
    txtfile << " Z -> mu mu" << endl;
    txtfile << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
    txtfile << "  pT > " << PT_CUT << endl;
    txtfile << "  |eta| < " << ETA_CUT << endl;
    txtfile << endl;

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        txtfile << "    *** Acceptance ***" << endl;
        txtfile << "          nominal: " << setw(12) << nSelv[ifile] << " / " << nEvtsv[ifile] << " = " << accv[ifile] << " +/- " << accErrv[ifile] << endl;
        txtfile << "  ==total efficiency==> " << setw(4) << accCorrv[ifile] / accv[ifile] << endl;
        txtfile << "     SF corrected: " << accCorrv[ifile] << " +/- " << accErrCorrv[ifile] << endl;
        txtfile << "          FSR unc: " << accCorrvFSR[ifile] << " / Sel: " << accCorrvFSR_I[ifile] << " / Sta: " << accCorrvFSR_S[ifile] << endl;
        txtfile << "           MC unc: " << accCorrvMC[ifile] << " / Sel: " << accCorrvMC_I[ifile] << " / Sta: " << accCorrvMC_S[ifile] << endl;
        txtfile << "          Bkg unc: " << accCorrvBkg[ifile] << " / Sel: " << accCorrvBkg_I[ifile] << " / Sta: " << accCorrvBkg_S[ifile] << endl;
        txtfile << "          Tag unc: " << accCorrvTag[ifile] << " / Sel: " << accCorrvTag_I[ifile] << " / Sta: " << accCorrvTag_S[ifile] << endl;
        // txtfile << "         Stat unc: " << accCorrvStat[ifile] << " / Sel: " << accCorrvStat_I[ifile]<< " / Sta: " << accCorrvStat_S[ifile] << endl;
        txtfile << "  fraction passing gen cut: " << nEvtsv[ifile] << " / " << nAllv[ifile] << " = " << nEvtsv[ifile] / nAllv[ifile] << endl;
        txtfile << endl;
    }
    txtfile.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sel_nums_only.txt", outputDir.Data());
    ofstream txtfile2;
    txtfile2.open(txtfname);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        txtfile2 << accCorrv[ifile] << " " << accErrCorrv[ifile] << endl;
        txtfile2 << accCorrvFSR[ifile] << " " << accCorrvFSR_I[ifile] << " " << accCorrvFSR_S[ifile] << endl;
        txtfile2 << accCorrvMC[ifile] << " " << accCorrvMC_I[ifile] << " " << accCorrvMC_S[ifile] << endl;
        txtfile2 << accCorrvBkg[ifile] << " " << accCorrvBkg_I[ifile] << " " << accCorrvBkg_S[ifile] << endl;
        txtfile2 << accCorrvTag[ifile] << " " << accCorrvTag_I[ifile] << " " << accCorrvTag_S[ifile] << endl;
        // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

        txtfile2 << endl;
    }
    txtfile2.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sit_unc.txt", outputDir.Data());
    ofstream txtfile3;
    txtfile3.open(txtfname);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        txtfile3 << accCorrv[ifile] << endl;
        txtfile3 << accCorrvFSR_I[ifile] << endl;
        txtfile3 << accCorrvMC_I[ifile] << endl;
        txtfile3 << accCorrvBkg_I[ifile] << endl;
        txtfile3 << accCorrvTag_I[ifile] << endl;
        // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

        txtfile3 << endl;
    }
    txtfile3.close();

    // char txtfname[100];
    sprintf(txtfname, "%s/sta_unc.txt", outputDir.Data());
    ofstream txtfile4;
    txtfile4.open(txtfname);

    //for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {
    for (UInt_t ifile = 0; ifile < 1; ifile++) {
        txtfile4 << accCorrv[ifile] << endl;
        txtfile4 << accCorrvFSR_S[ifile] << endl;
        txtfile4 << accCorrvMC_S[ifile] << endl;
        txtfile4 << accCorrvBkg_S[ifile] << endl;
        txtfile4 << accCorrvTag_S[ifile] << endl;
        // txtfile << accCorrvFSR[ifile]  << ", " << accCorrvMC[ifile] << ", " << accCorrvBkg[ifile] << ", " << accCorrvTag[ifile] << endl;

        txtfile4 << endl;
    }
    txtfile4.close();

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;

    gBenchmark->Show("computeAccSelZmm");
}
