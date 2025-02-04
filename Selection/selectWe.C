//================================================================================================
//
// Select W->enu candidates
//
//  * outputs ROOT files of events passing selection
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLorentzVector.h" // 4-vector class
#include "TRandom.h"
#include <TBenchmark.h>   // class to track macro running statistics
#include <TClonesArray.h> // ROOT array class
#include <TFile.h>        // file handle class
#include <TMath.h>        // ROOT math library
#include <TROOT.h>        // access to gROOT, entry point to ROOT system
#include <TSystem.h>      // interface to OS
#include <TTree.h>        // class to access ntuples
#include <TVector2.h>     // 2D vector class
#include <fstream>        // functions for file I/O
#include <iomanip>        // functions to format standard I/O
#include <iostream>       // standard I/O
#include <vector>         // STL vector class

// define structures to read in ntuple

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

// lumi section selection with JSON files
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"

#include "MitEwk13TeV/EleScale/EnergyScaleCorrection.h" //EGMSmear
#include "MitEwk13TeV/Selection/ConfParse.hh"           // input conf file parser
#include "MitEwk13TeV/Utils/CSample.hh"                 // helper class to handle samples
#include "MitEwk13TeV/Utils/LeptonCorr.hh"              // electron scale and resolution corrections
#include "MitEwk13TeV/Utils/LeptonIDCuts.hh"            // helper functions for lepton ID selection
#include "MitEwk13TeV/Utils/MyTools.hh"                 // various helper functions
#include "MitEwk13TeV/Utils/PrefiringEfficiency.cc"     // prefiring efficiency functions
#endif

//=== MAIN MACRO =================================================================================================

void selectWe(const TString conf = "we.conf", // input file
              const TString outputDir = ".",  // output directory
              const Bool_t doScaleCorr = 0,   // apply energy scale corrections?
              const Int_t sigma = 0,
              const Bool_t doPU = 0,
              const Bool_t is13TeV = 1,
              const Int_t NSEC = 1,
              const Int_t ITH = 0,
              const Int_t ISample = -1)
{
    gBenchmark->Start("selectWe");

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================
    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;

    const Double_t VETO_PT = 10;
    const Double_t VETO_ETA = 2.4;

    const Double_t ECAL_GAP_LOW = 1.4442;
    const Double_t ECAL_GAP_HIGH = 1.566;

    const Double_t ELE_MASS = 0.000511;
    const Int_t BOSON_ID = 24;
    const Int_t LEPTON_ID = 11;

    // theory unc variations, taken from
    // https://github.com/yongbinfeng/BaconProd/blob/master/Ntupler/src/FillerGenInfo.cc#L45-L48
    // 8 qcd scale + 1 default pdf +100 pdf + 2 alphaS
    // for W's the weights are a bit different. there are
    // one extra redundant qcd scale variation (always 1)
    // so it becomes 9 qcd scale + 1 default pdf + 100 pdf + 1alphaS
    // the other alphaS is not saved in the ntuples.
    const Int_t nTHEORYUNC = 111;

    const TString envStr = (TString)gSystem->Getenv("CMSSW_BASE") + "/src/";

    // load trigger menu
    const baconhep::TTrigger triggerMenu((envStr + "BaconAna/DataFormats/data/HLT_50nsGRun").Data());

    const TString prefireEcalFileName = envStr + "MitEwk13TeV/Utils/All2017Gand2017HPrefiringMaps.root";
    const TString prefireMuonFileName = envStr + "MitEwk13TeV/Utils/L1MuonPrefiringParametriations.root";
    PrefiringEfficiency pfire(prefireEcalFileName.Data(), (is13TeV ? "2017H" : "2017G"), prefireMuonFileName.Data());

    // load pileup reweighting file
    TFile *f_rw = TFile::Open(envStr + "MitEwk13TeV/Tools/puWeights_76x.root", "read");
    TH1D *h_rw = (TH1D *)f_rw->Get("puWeights");
    TH1D *h_rw_up = (TH1D *)f_rw->Get("puWeightsUp");
    TH1D *h_rw_down = (TH1D *)f_rw->Get("puWeightsDown");

    const TString corrFiles = envStr + "MitEwk13TeV/EleScale/Run2017_LowPU_v2";
    EnergyScaleCorrection eleCorr(corrFiles.Data(), EnergyScaleCorrection::ECALELF);

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    vector<TString> snamev;    // sample name (for output files)
    vector<CSample *> samplev; // data/MC samples

    // parse .conf file
    confParse(conf, snamev, samplev);
    const Bool_t hasData = (samplev[0]->fnamev.size() > 0);

    // Create output directory
    gSystem->mkdir(outputDir, kTRUE);
    const TString ntupDir = outputDir + TString("/ntuples");
    // const TString ntupDir = outputDir + TString("/ntuples_") + Form("%d", ITH) + TString("_") + Form("%d", NSEC);
    gSystem->mkdir(ntupDir, kTRUE);

    // Declare output ntuple variables
    UInt_t runNum, lumiSec, evtNum;
    UInt_t npv, npu;
    TLorentzVector *genV = 0, *genLep = 0, *genNu = 0;
    Float_t scale1fb, scale1fbUp, scale1fbDown;
    // these are naming conventions
    // actually the prefireWeight here should be the non-prefirable probability
    Float_t prefireEcal = 1, prefireEcalUp = 1, prefireEcalDown = 1;
    Float_t prefirePhoton = 1, prefirePhotUp = 1, prefirePhotDown = 1;
    Float_t prefireJet = 1, prefireJetUp = 1, prefireJetDown = 1;
    Float_t prefireMuon = 1, prefireMuonUp = 1, prefireMuonDown = 1, prefireMuonStatUp = 1, prefireMuonStatDown = 1, prefireMuonSystUp = 1, prefireMuonSystDown = 1;
    Float_t prefireWeight = 1, prefireUp = 1, prefireDown = 1;
    Float_t met, metPhi;           //, mt, u1, u2;
    Float_t puppiMet, puppiMetPhi; //, puppiMt, puppiU1, puppiU2;
    Int_t q;
    TLorentzVector *lep = 0, *lep_raw = 0;
    Float_t lepError = 0;
    Float_t pfChIso, pfGamIso, pfNeuIso;
    Float_t pfCombIso;
    TLorentzVector *sc = 0;
    vector<Double_t> lheweight(nTHEORYUNC, 1.0);

    // Data structures to store info from TTrees
    baconhep::TEventInfo *info = new baconhep::TEventInfo();
    baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
    TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
    TClonesArray *electronArr = new TClonesArray("baconhep::TElectron");
    TClonesArray *muonArr = new TClonesArray("baconhep::TMuon");
    TClonesArray *scArr = new TClonesArray("baconhep::TPhoton");
    TClonesArray *vertexArr = new TClonesArray("baconhep::TVertex");
    TClonesArray *jetArr = new TClonesArray("baconhep::TJet");

    TFile *infile = 0;
    TTree *eventTree = 0;

    // loop over samples
    for (UInt_t isam = 0; isam < samplev.size(); isam++)
    {
        // std::cout << "ISample " << ISample << std::endl;
        if ((ISample >= 0) && ((int)isam != ISample))
            // only run the ISample in the config
            continue;
        std::cout << "start running for isamp " << isam << std::endl;

        // Assume data sample is first sample in .conf file
        // If sample is empty (i.e. contains no ntuple files), skip to next sample
        Bool_t isData = false;
        // if (isam == 0 && !hasData)
        //     continue;
        // else if (isam == 0)
        //     isData = kTRUE;
        if (snamev[isam].CompareTo("data", TString::kIgnoreCase) == 0)
            isData = true;

        Bool_t isSignal = (snamev[isam].Contains("we"));
        Bool_t isWrongFlavor = (snamev[isam].Contains("wx"));

        Bool_t isRecoil = (snamev[isam].Contains("zxx") || isSignal || isWrongFlavor);
        Bool_t noGen = (snamev[isam].Contains("zz") || snamev[isam].Contains("wz") || snamev[isam].Contains("ww") || snamev[isam].Contains("zz2l") || snamev[isam].Contains("zz4l"));

        Bool_t isDY = (snamev[isam].CompareTo("zxx", TString::kIgnoreCase) == 0);

        CSample *samp = samplev[isam];

        // Set up output ntuple
        TString outfilename = ntupDir + TString("/") + snamev[isam];
        if (isam != 0 && !doScaleCorr)
            // not data and no scale correction
            outfilename += TString("_select.raw");
        else
            outfilename += TString("_select");
        if (NSEC != 1)
        {
            outfilename += TString("_") + Form("%dSections", NSEC) + TString("_") + Form("%d", ITH);
        }
        outfilename += TString(".root");
        cout << outfilename << endl;

        for (unsigned itheory = 0; itheory < nTHEORYUNC; itheory++)
        {
            lheweight[itheory] = 1.0;
        }

        TFile *outFile = new TFile(outfilename, "RECREATE");
        TTree *outTree = new TTree("Events", "Events");
        outTree->Branch("runNum", &runNum, "runNum/i");       // event run number
        outTree->Branch("lumiSec", &lumiSec, "lumiSec/i");    // event lumi section
        outTree->Branch("evtNum", &evtNum, "evtNum/i");       // event number
        outTree->Branch("npv", &npv, "npv/i");                // number of primary vertices
        outTree->Branch("npu", &npu, "npu/i");                // number of in-time PU events (MC)
        outTree->Branch("genV", "TLorentzVector", &genV);     // GEN boson 4-vector (signal MC)
        outTree->Branch("genLep", "TLorentzVector", &genLep); // GEN lepton 4-vector (signal MC)
        outTree->Branch("genNu", "TLorentzVector", &genNu);   // GEN lepton 4-vector (signal MC)
        outTree->Branch("prefireEcal", &prefireEcal, "prefireEcal/F");
        outTree->Branch("prefireEcalUp", &prefireEcalUp, "prefireEcalUp/F");
        outTree->Branch("prefireEcalDown", &prefireEcalDown, "prefireEcalDown/F");
        outTree->Branch("prefirePhoton", &prefirePhoton, "prefirePhoton/F");
        outTree->Branch("prefirePhotUp", &prefirePhotUp, "prefirePhotUp/F");
        outTree->Branch("prefirePhotDown", &prefirePhotDown, "prefirePhotDown/F");
        outTree->Branch("prefireJet", &prefireJet, "prefireJet/F");
        outTree->Branch("prefireJetUp", &prefireJetUp, "prefireJetUp/F");
        outTree->Branch("prefireJetDown", &prefireJetDown, "prefireJetDown/F");
        outTree->Branch("prefireMuon", &prefireMuon, "prefireMuon/F");
        outTree->Branch("prefireMuonUp", &prefireMuonUp, "prefireMuonUp/F");
        outTree->Branch("prefireMuonDown", &prefireMuonDown, "prefireMuonDown/F");
        outTree->Branch("prefireMuonStatUp", &prefireMuonStatUp, "prefireMuonStatUp/F");
        outTree->Branch("prefireMuonStatDown", &prefireMuonStatDown, "prefireMuonStatDown/F");
        outTree->Branch("prefireMuonSystUp", &prefireMuonSystUp, "prefireMuonSystUp/F");
        outTree->Branch("prefireMuonSystDown", &prefireMuonSystDown, "prefireMuonSystDown/F");
        outTree->Branch("prefireWeight", &prefireWeight, "prefireWeight/F");
        outTree->Branch("prefireUp", &prefireUp, "prefireUp/F");
        outTree->Branch("prefireDown", &prefireDown, "prefireDown/F");
        outTree->Branch("scale1fb", &scale1fb, "scale1fb/F");             // event weight per 1/fb (MC)
        outTree->Branch("scale1fbUp", &scale1fbUp, "scale1fbUp/F");       // event weight per 1/fb (MC)
        outTree->Branch("scale1fbDown", &scale1fbDown, "scale1fbDown/F"); // event weight per 1/fb (MC)
        outTree->Branch("met", &met, "met/F");                            // MET
        outTree->Branch("metPhi", &metPhi, "metPhi/F");                   // phi(MET)
        outTree->Branch("puppiMet", &puppiMet, "puppiMet/F");             // Puppi MET
        outTree->Branch("puppiMetPhi", &puppiMetPhi, "puppiMetPhi/F");    // phi(Puppi MET)
        outTree->Branch("q", &q, "q/I");                                  // lepton charge
        outTree->Branch("lep", "TLorentzVector", &lep);                   // lepton 4-vector
        outTree->Branch("lep_raw", "TLorentzVector", &lep_raw);           // lepton 4-vector
        outTree->Branch("lepError", &lepError, "lepError/F");             // track isolation of tag lepton
        outTree->Branch("pfChIso", &pfChIso, "pfChIso/F");                // PF charged hadron isolation of lepton
        outTree->Branch("pfGamIso", &pfGamIso, "pfGamIso/F");             // PF photon isolation of lepton
        outTree->Branch("pfNeuIso", &pfNeuIso, "pfNeuIso/F");             // PF neutral hadron isolation of lepton
        outTree->Branch("pfCombIso", &pfCombIso, "pfCombIso/F");          // PF combined isolation of electron
        outTree->Branch("sc", "TLorentzVector", &sc);                     // supercluster 4-vector
        outTree->Branch("lheweight", "vector<Double_t>", &lheweight);     // LHE weights

        TH1D *hGenWeights = new TH1D("hGenWeights", "hGenWeights", 2, -1., 1.);
        // save the sum of weights with different LHE weight variations (QCD scale, pdf, alphaS)
        TH1D *hLHEWeightSum = new TH1D("hLHEWeightSum", "hLHEWeightSum", nTHEORYUNC, 0, nTHEORYUNC);

        //
        // loop through files
        //
        const UInt_t nfiles = samp->fnamev.size();
        for (UInt_t ifile = 0; ifile < nfiles; ifile++)
        {
            // Read input file and get the TTrees
            cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... ";
            cout.flush();
            infile = TFile::Open(samp->fnamev[ifile]);
            assert(infile);

            Bool_t hasJSON = kFALSE;
            baconhep::RunLumiRangeMap rlrm;
            if (!samp->jsonv[ifile].Contains("NONE"))
            {
                hasJSON = kTRUE;
                rlrm.addJSONFile(samp->jsonv[ifile].Data());
            }

            eventTree = (TTree *)infile->Get("Events");
            assert(eventTree);
            Bool_t hasJet = eventTree->GetBranchStatus("AK4");
            eventTree->SetBranchAddress("Info", &info);
            TBranch *infoBr = eventTree->GetBranch("Info");
            eventTree->SetBranchAddress("Electron", &electronArr);
            TBranch *electronBr = eventTree->GetBranch("Electron");
            eventTree->SetBranchAddress("Muon", &muonArr);
            TBranch *muonBr = eventTree->GetBranch("Muon");
            eventTree->SetBranchAddress("PV", &vertexArr);
            TBranch *vertexBr = eventTree->GetBranch("PV");
            eventTree->SetBranchAddress("Photon", &scArr);
            TBranch *scBr = eventTree->GetBranch("Photon");
            if (hasJet)
                eventTree->SetBranchAddress("AK4", &jetArr);
            TBranch *jetBr = eventTree->GetBranch("AK4");

            Bool_t hasGen = (eventTree->GetBranchStatus("GenEvtInfo") && !noGen);
            TBranch *genBr = 0, *genPartBr = 0;
            if (hasGen)
            {
                eventTree->SetBranchAddress("GenEvtInfo", &gen);
                genBr = eventTree->GetBranch("GenEvtInfo");
                eventTree->SetBranchAddress("GenParticle", &genPartArr);
                genPartBr = eventTree->GetBranch("GenParticle");
            }
            // Compute MC event weight per 1/fb
            const Double_t xsec = samp->xsecv[ifile];
            Double_t puWeight = 0, puWeightUp = 0, puWeightDown = 0;

            //////////////////////////////////////////////////////////
            // Real selection loop
            //
            // loop over events
            //
            Long64_t nevents = eventTree->GetEntries();
            Long64_t IBEGIN = 0;
            Long64_t IEND = nevents;

            if (NSEC != 1)
            {
                cout << "n sections " << NSEC << " ith part " << ITH << endl;
                Long64_t nevents_per_job = nevents / NSEC;
                Long64_t remainder = nevents % NSEC;
                IBEGIN = ITH * nevents_per_job;
                IEND = (ITH + 1) * nevents_per_job;
                if (ITH == NSEC - 1)
                {
                    // for the last section, add the remaining into these;
                    IEND = nevents;
                }
                cout << "start, end events: " << IBEGIN << " " << IEND << endl;
            }

            std::cout << "Number of Events = " << eventTree->GetEntries() << ", among these processing " << IEND - IBEGIN << std::endl;

            Double_t nsel = 0, nselvar = 0;
            for (UInt_t ientry = IBEGIN; ientry < IEND; ientry++)
            {
                infoBr->GetEntry(ientry);

                int printIndex = (int)(eventTree->GetEntries() * 0.01);
                if (ientry % printIndex == 0)
                    cout << "Processing event " << ientry << ". " << (int)(100 * (ientry / (double)eventTree->GetEntries())) << " percent done with this file." << endl;

                Double_t weight = xsec, weightUp = xsec, weightDown = xsec;
                if (hasGen)
                {
                    genPartArr->Clear();
                    genBr->GetEntry(ientry);
                    genPartBr->GetEntry(ientry);
                    puWeight = doPU ? h_rw->GetBinContent(h_rw->FindBin(info->nPUmean)) : 1.;
                    puWeightUp = doPU ? h_rw_up->GetBinContent(h_rw_up->FindBin(info->nPUmean)) : 1.;
                    puWeightDown = doPU ? h_rw_down->GetBinContent(h_rw_down->FindBin(info->nPUmean)) : 1.;
                    int genweight = gen->weight > 0 ? 1 : -1;
                    hGenWeights->Fill(0.0, genweight);
                    weight *= genweight * puWeight;
                    weightUp *= genweight * puWeightUp;
                    weightDown *= genweight * puWeightDown;
                    for (unsigned itheory = 0; itheory < nTHEORYUNC; itheory++)
                    {
                        lheweight[itheory] = gen->lheweight[itheory];
                        hLHEWeightSum->Fill(itheory, genweight * gen->lheweight[itheory]);
                    }
                }
                else
                {
                    hGenWeights->Fill(0.0, 1.0);
                }
                // veto w -> xv decays for signal and w -> mv for bacground samples (needed for inclusive WToLNu sample)
                if (isWrongFlavor && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID)) == LEPTON_ID)
                    continue;
                else if (isSignal && hasGen && fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                    continue;
                // check for certified lumi (if applicable)
                baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
                if (hasJSON && !rlrm.hasRunLumi(rl))
                    continue;
                // trigger requirement
                if (!isEleTrigger(triggerMenu, info->triggerBits, isData, is13TeV))
                    continue;
                // good vertex requirement
                if (!(info->hasGoodPV))
                    continue;
                //
                // SELECTION PROCEDURE:
                //  (1) Look for 1 good electron matched to trigger
                //  (2) Reject event if another electron is present passing looser cuts
                //
                electronArr->Clear();
                electronBr->GetEntry(ientry);
                scArr->Clear();
                scBr->GetEntry(ientry);
                muonArr->Clear();
                muonBr->GetEntry(ientry);
                jetArr->Clear();
                if (hasJet)
                    jetBr->GetEntry(ientry);

                Int_t nLooseLep = 0;
                const baconhep::TElectron *goodEle = 0;
                TLorentzVector vEle(0, 0, 0, 0);
                TLorentzVector vGoodEle(0, 0, 0, 0);
                Bool_t passSel = kFALSE;
                double eleRamdom = gRandom->Gaus(0, 1);

                for (Int_t i = 0; i < muonArr->GetEntriesFast(); i++)
                {
                    const baconhep::TMuon *mu = (baconhep::TMuon *)((*muonArr)[i]);
                    if (fabs(mu->eta) > VETO_ETA)
                        continue; // loose lepton |eta| cut
                    if (mu->pt < VETO_PT)
                        continue; // loose lepton pT cut
                    if (passMuonLooseID(mu))
                        nLooseLep++; // loose lepton selection
                    if (nLooseLep > 0)
                    { // extra lepton veto
                        passSel = kFALSE;
                        break;
                    }
                }

                for (Int_t i = 0; i < electronArr->GetEntriesFast(); i++)
                {
                    const baconhep::TElectron *ele = (baconhep::TElectron *)((*electronArr)[i]);
                    vEle.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, ELE_MASS);
                    // check ECAL gap
                    // if(fabs(vEle.Eta())>=ECAL_GAP_LOW && fabs(vEle.Eta())<=ECAL_GAP_HIGH) continue;

                    if (doScaleCorr && (ele->r9 < 1.))
                    {
                        float eleAbsEta = fabs(vEle.Eta());
                        double eTregress = ele->ecalEnergy / cosh(fabs(ele->eta));
                        if (snamev[isam].Contains("data"))
                        { // Data
                            int runNumber = is13TeV ? info->runNum : 306936;
                            float eleScale = eleCorr.scaleCorr(runNumber, eTregress, eleAbsEta, ele->r9);
                            float eleError = eleCorr.scaleCorrUncert(runNumber, eTregress, eleAbsEta, ele->r9);
                            (vEle) *= eleScale * (1 + sigma * eleError);
                            lepError = eleError;
                        }
                        else
                        { // MC
                            float eleSmear = eleCorr.smearingSigma(info->runNum, eTregress, eleAbsEta, ele->r9, 12, sigma, 0.);
                            (vEle) *= 1. + eleSmear * eleRamdom;
                            float eleSmearEP = eleCorr.smearingSigma(info->runNum, eTregress, eleAbsEta, ele->r9, 12, 1., 0.);
                            float eleSmearEM = eleCorr.smearingSigma(info->runNum, eTregress, eleAbsEta, ele->r9, 12, -1., 0.);
                            lepError = eleRamdom * std::hypot(eleSmearEP - eleSmear, eleSmearEM - eleSmear);
                        }
                    }

                    // apply scale and resolution corrections to MC
                    if (fabs(vEle.Eta()) > VETO_ETA)
                        continue;
                    if (vEle.Pt() < VETO_PT)
                        continue;
                    if (fabs(vEle.Eta()) >= ECAL_GAP_LOW && fabs(vEle.Eta()) <= ECAL_GAP_HIGH)
                        continue;
                    if (passEleLooseID(ele, vEle, info->rhoIso))
                        nLooseLep++;
                    if (nLooseLep > 1)
                    { // extra lepton veto
                        passSel = kFALSE;
                        break;
                    }
                    if (vEle.Pt() < PT_CUT)
                        continue; // lepton pT cut
                    if (fabs(vEle.Eta()) > ETA_CUT)
                        continue; // lepton |eta| cut
                    if (!passEleMediumID(ele, vEle, info->rhoIso))
                        continue; // lepton selection
                    if (!isEleTriggerObj(triggerMenu, ele->hltMatchBits, kFALSE, isData, is13TeV))
                        continue;
                    passSel = kTRUE;
                    goodEle = ele;
                    vGoodEle = vEle;
                }

                if (passSel)
                {
                    //******* We have a W candidate! HURRAY! ********
                    nsel += isData ? 1 : weight;
                    nselvar += isData ? 1 : weight * weight;

                    if (!isData)
                    {
                        pfire.setObjects(scArr, jetArr, muonArr);
                        pfire.computePhotonsOnly(prefirePhoton, prefirePhotUp, prefirePhotDown);
                        pfire.computeJetsOnly(prefireJet, prefireJetUp, prefireJetDown);
                        pfire.computeEcalsOnly(prefireEcal, prefireEcalUp, prefireEcalDown);
                        pfire.computeMuonsOnly(prefireMuon, prefireMuonUp, prefireMuonDown, prefireMuonStatUp, prefireMuonStatDown, prefireMuonSystUp, prefireMuonSystDown);
                        prefireWeight = prefireEcal * prefireMuon;
                        prefireUp = prefireEcalUp * prefireMuonUp;
                        prefireDown = prefireEcalDown * prefireMuonDown;
                    }

                    TLorentzVector vLep(0, 0, 0, 0);
                    TLorentzVector vSC(0, 0, 0, 0);
                    TLorentzVector vLep_raw(0, 0, 0, 0);
                    vLep = vGoodEle;
                    vLep_raw.SetPtEtaPhiM(goodEle->pt, goodEle->eta, goodEle->phi, ELE_MASS);

                    //
                    // Fill tree
                    //
                    runNum = info->runNum;
                    lumiSec = info->lumiSec;
                    evtNum = info->evtNum;

                    vertexArr->Clear();
                    vertexBr->GetEntry(ientry);

                    npv = vertexArr->GetEntries();
                    npu = info->nPUmean;
                    genV = new TLorentzVector(0, 0, 0, 0);
                    genLep = new TLorentzVector(0, 0, 0, 0);
                    genNu = new TLorentzVector(0, 0, 0, 0);

                    if (isRecoil && hasGen)
                    {
                        TLorentzVector *gvec = new TLorentzVector(0, 0, 0, 0);
                        TLorentzVector *glep1 = new TLorentzVector(0, 0, 0, 0);
                        TLorentzVector *glep2 = new TLorentzVector(0, 0, 0, 0);
                        TLorentzVector *glepB1 = new TLorentzVector(0, 0, 0, 0);
                        TLorentzVector *glepB2 = new TLorentzVector(0, 0, 0, 0);
                        if (isDY)
                            toolbox::fillGenBorn(genPartArr, 23, gvec, glepB1, glepB2, glep1, glep2);
                        else
                            toolbox::fillGenBorn(genPartArr, BOSON_ID, gvec, glepB1, glepB2, glep1, glep2);

                        // TLorentzVector tvec = *glep1 + *glep2;
                        genV = new TLorentzVector(0, 0, 0, 0);
                        // genV->SetPtEtaPhiM(tvec.Pt(), tvec.Eta(), tvec.Phi(), tvec.M());
                        genV->SetPtEtaPhiM(gvec->Pt(), gvec->Eta(), gvec->Phi(), gvec->M());
                        if (gvec && glep1)
                        {
                            if (!isDY)
                            {
                                if (toolbox::flavor(genPartArr, BOSON_ID) < 0)
                                { // means it's a W+ and charged lepton is anti-particle
                                    genLep->SetPtEtaPhiM(glep2->Pt(), glep2->Eta(), glep2->Phi(), glep2->M());
                                    genNu->SetPtEtaPhiM(glep1->Pt(), glep1->Eta(), glep1->Phi(), glep1->M());
                                }
                                if (toolbox::flavor(genPartArr, BOSON_ID) > 0)
                                { // means it's a W- and charged lepton is particle
                                    genLep->SetPtEtaPhiM(glep1->Pt(), glep1->Eta(), glep1->Phi(), glep1->M());
                                    genNu->SetPtEtaPhiM(glep2->Pt(), glep2->Eta(), glep2->Phi(), glep2->M());
                                }
                            }
                            else
                            {
                                if (goodEle->q > 0)
                                {
                                    // lepton pdgId: -11
                                    genLep->SetPtEtaPhiM(glep2->Pt(), glep2->Eta(), glep2->Phi(), glep2->M());
                                    genNu->SetPtEtaPhiM(glep1->Pt(), glep1->Eta(), glep1->Phi(), glep1->M());
                                }
                                else
                                {
                                    // lepton pdgId: 11
                                    genLep->SetPtEtaPhiM(glep1->Pt(), glep1->Eta(), glep1->Phi(), glep1->M());
                                    genNu->SetPtEtaPhiM(glep2->Pt(), glep2->Eta(), glep2->Phi(), glep2->M());
                                }
                            }
                        }
                        delete gvec;
                        delete glep1;
                        delete glep2;
                        delete glepB1;
                        delete glepB2;
                        gvec = 0;
                        glep1 = 0;
                        glep2 = 0;
                        glepB1 = 0;
                        glepB2 = 0;
                    }
                    scale1fb = weight;
                    scale1fbUp = weightUp;
                    scale1fbDown = weightDown;

                    met = info->pfMETC;
                    metPhi = info->pfMETCphi;
                    puppiMet = info->puppET;
                    puppiMetPhi = info->puppETphi;
                    q = goodEle->q;
                    lep = &vLep;
                    lep_raw = &vLep_raw;

                    ///// electron specific /////
                    sc = &vSC;
                    pfChIso = goodEle->chHadIso;
                    pfGamIso = goodEle->gammaIso;
                    pfNeuIso = goodEle->neuHadIso;
                    pfCombIso = goodEle->chHadIso + TMath::Max(goodEle->neuHadIso + goodEle->gammaIso - (info->rhoIso) * getEffAreaEl(goodEle->eta), 0.);

                    outTree->Fill();
                    delete genV;
                    delete genLep;
                    delete genNu;
                    genV = 0, genLep = 0, lep = 0, lep_raw = 0, genNu = 0, sc = 0;
                    // reset everything to 1
                    prefirePhoton = 1;
                    prefirePhotUp = 1;
                    prefirePhotDown = 1;
                    prefireJet = 1;
                    prefireJetUp = 1;
                    prefireJetDown = 1;
                    prefireEcal = 1;
                    prefireEcalUp = 1;
                    prefireEcalDown = 1;
                    prefireMuon = 1;
                    prefireMuonUp = 1;
                    prefireMuonDown = 1;
                    prefireMuonStatUp = 1;
                    prefireMuonStatDown = 1;
                    prefireMuonSystUp = 1;
                    prefireMuonSystDown = 1;
                    prefireWeight = 1;
                    prefireUp = 1;
                    prefireDown = 1;
                }
            }
            delete infile;
            infile = 0, eventTree = 0;

            cout << nsel << " +/- " << sqrt(nselvar);
            if (isam != 0)
                cout << " per 1/pb";
            cout << endl;
        }
        outFile->cd();
        hGenWeights->Write();
        hLHEWeightSum->Write();
        outFile->Write();
        outFile->Close();
    }
    delete h_rw;
    delete h_rw_up;
    delete h_rw_down;
    delete f_rw;
    delete info;
    delete gen;
    delete genPartArr;
    delete electronArr;
    delete vertexArr;

    //--------------------------------------------------------------------------------------------------------------
    // Output
    //==============================================================================================================

    cout << "*" << endl;
    cout << "* SUMMARY" << endl;
    cout << "*--------------------------------------------------" << endl;
    cout << " W -> e nu" << endl;
    cout << "  pT > " << PT_CUT << endl;
    cout << "  |eta| < " << ETA_CUT << endl;
    if (doScaleCorr)
        cout << "  *** Scale corrections applied ***" << endl;
    cout << endl;

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;

    gBenchmark->Show("selectWe");
}
