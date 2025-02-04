//================================================================================================
//
//  Perform recoil corrections, rochester corrections, make new branches for the info
//
//  * outputs another ntuple but with new branches
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TLorentzVector.h" // 4-vector class
#include <TBenchmark.h>     // class to track macro running statistics
#include <TFile.h>          // file handle class
#include <TGaxis.h>
#include <TH1D.h>  // histogram class
#include <TROOT.h> // access to gROOT, entry point to ROOT system
#include <TRandom3.h>
#include <TStyle.h>  // class to handle ROOT plotting styles
#include <TSystem.h> // interface to OS
#include <TTree.h>   // class to access ntuples
#include <fstream>   // functions for file I/O
#include <iomanip>   // functions to format standard I/O
#include <iostream>  // standard I/O
#include <sstream>   // class for parsing strings
#include <string>    // C++ string class
#include <vector>    // STL vector class

#include "MitEwk13TeV/Utils/MyTools.hh" // various helper functions
#include "MitEwk13TeV/Utils/RecoilCorrector.hh"
#include "MitEwk13TeV/Utils/METXYCorrector.hh"
#include "MitEwk13TeV/Utils/AppEffSF.cc"

#endif

//=== MAIN MACRO =================================================================================================

void eleNtupleMod(const TString outputDir, // output directory
                  const TString inputDir,  // input directory
                  const TString sqrts,     // 13 or 5 TeV string specifier
                  const TString fileName,  // both the input and output final file name i.e. data_select.root
                  const Int_t NSEC = 1,
                  const Int_t ITH = 0)
{
    std::cout << "start running new " << std::endl;
    gBenchmark->Start("fitWe");
    gROOT->SetBatch(1);

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================
    // flage to control applying the recoil corrections
    bool doInclusive = true; // This should be the standard recoil correction: 3-gaussian inclusive eta
    bool doKeys = true;      // RooKeysPDF instead of 3-Gaus
    bool doEta = true;       // eta-binned 3-Gaus fit
    bool doStat = false;     //  Statistical Uncertainty
    int nNV = 10;
    bool doPF = true;

    // std::string u1_name;
    // std::string u2_name;
    std::string met_name;
    std::string metPhi_name;
    //   std::string recoilType;

    // Control the types of uncertainties
    enum
    {
        no,
        cent,
        eta,
        keys,
        ru,
        rd,
        stat0,
        stat1,
        stat2,
        stat3,
        stat4,
        stat5,
        stat6,
        stat7,
        stat8,
        stat9,
        org
    };
    const string vMET[] = {"no", "cent", "eta", "keys", "ru", "rd", "stat0", "stat1", "stat2", "stat3", "stat4", "stat5", "stat6", "stat7", "stat8", "stat9", "org"};
    int nMET = sizeof(vMET) / sizeof(vMET[0]);
    // int ns = nMET - nNV;
    int ns = 6;
    // front half should be nMET-nNV

    enum
    {
        main,
        mc,
        fsr,
        bkg,
        tagpt,
        effstat,
        pfireu,
        pfired,
        pfireecalu,
        pfireecald,
        pfiremuu,
        pfiremud,
        effstat_lepPos,
        effstat_lepNeg
    };
    const string vWeight[] = {"eff", "mc", "fsr", "bkg", "tagpt", "effstat", "pfireu", "pfired", "pfireecalu", "pfireecald", "pfiremuu", "pfiremud", "effstat_lepPos", "effstat_lepNeg"};
    int nWeight = sizeof(vWeight) / sizeof(vWeight[0]);

    // u1_name = "u1";
    // u2_name = "u2";
    met_name = "met";
    metPhi_name = "metPhi";

    // don't think these are really necessary but leaving them for now

    const Double_t ECAL_GAP_LOW = 1.4442;
    const Double_t ECAL_GAP_HIGH = 1.566;
    // const Double_t ECAL_GAP_LOW  = 10.;
    // const Double_t ECAL_GAP_HIGH = 10.;

    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;
    const Double_t ELE_MASS = 0.000511;

    const bool doMETXYCorrection = true;

    const TString envStr = (TString)gSystem->Getenv("CMSSW_BASE") + "/src/";

    TString baseDir = envStr + "Corrections/Efficiency/" + sqrts + "/results/Zee/";
    AppEffSF effs(baseDir);
    effs.loadHLT("EleHLTEff_aMCxPythia", "Positive", "Negative");
    effs.loadSel("EleGSFSelEff_aMCxPythia", "Combined", "Combined");
    //
    // Warning: this needs to be updated for 5TeV
    //
    TString SysFileGSFSel = envStr + "Corrections/Efficiency/" + sqrts + "/results/Systematics/SysUnc_EleGSFSelEff.root";
    effs.loadUncSel(SysFileGSFSel);
    TH2D *hErr = new TH2D("hErr", "", 10, 0, 10, 20, 0, 20);

    Bool_t isData = (fileName.Contains("data_select"));
    std::cout << "isData ? " << isData << std::endl;

    Bool_t isRecoil = (fileName.Contains("we_select") || fileName.Contains("we0_select") || fileName.Contains("we1_select") || fileName.Contains("we2_select") || fileName.Contains("wx_select") || fileName.Contains("wx0_select") || fileName.Contains("wx1_select") || fileName.Contains("wx2_select") || fileName.Contains("zxx_select"));
    std::cout << "do Recoil " << isRecoil << std::endl;

    if (inputDir.Contains("Anti") && isRecoil)
    {
        // Anti isolated region is for QCD bkg estimations
        // no need to do recoil uncertainties
        doInclusive = true;
        doKeys = false;
        doEta = false;
        doStat = false;
    }

    if (isData || (!isRecoil))
    {
        doInclusive = false;
        doKeys = false;
        doEta = false;
        doStat = false;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------
    //   Load the Recoil Correction Files
    // ------------------------------------------------------------------------------------------------------------------------------------------
    // ===================== Recoil correction files ============================
    const TString directory((envStr + "Corrections/RecoilNew").Data());

    // New Recoil Correctors for everything
    RecoilCorrector *rcMainWp = new RecoilCorrector("", "");
    RecoilCorrector *rcMainWm = new RecoilCorrector("", "");
    vector<RecoilCorrector *> rcStatW;
    for (int i = 0; i < nNV; ++i)
    {
        RecoilCorrector *tempStatW = new RecoilCorrector("", "");
        rcStatW.push_back(tempStatW);
    }
    // RecoilCorrector *rcStatW     = new  RecoilCorrector("","");
    RecoilCorrector *rcKeysWp = new RecoilCorrector("", "");
    RecoilCorrector *rcKeysWm = new RecoilCorrector("", "");
    RecoilCorrector *rcEta05Wp = new RecoilCorrector("", "");
    RecoilCorrector *rcEta05Wm = new RecoilCorrector("", "");
    RecoilCorrector *rcEta051Wp = new RecoilCorrector("", "");
    RecoilCorrector *rcEta051Wm = new RecoilCorrector("", "");
    RecoilCorrector *rcEta1Wp = new RecoilCorrector("", "");
    RecoilCorrector *rcEta1Wm = new RecoilCorrector("", "");
    // also make sure to add the Wm stuff
    if (doInclusive)
    {
        rcMainWp->loadRooWorkspacesMCtoCorrect(Form("%s/WepMC_PF_%s_2G/", directory.Data(), sqrts.Data()));
        rcMainWp->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_2G/", directory.Data(), sqrts.Data()));
        rcMainWp->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_2G/", directory.Data(), sqrts.Data()));

        rcMainWm->loadRooWorkspacesMCtoCorrect(Form("%s/WemMC_PF_%s_2G/", directory.Data(), sqrts.Data()));
        rcMainWm->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_2G/", directory.Data(), sqrts.Data()));
        rcMainWm->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_2G/", directory.Data(), sqrts.Data()));
    }
    if (doStat)
    {
        int rec_sig = 1;
        for (int i = 0; i < nNV; i++)
        {
            rcStatW[i]->loadRooWorkspacesDiagMCtoCorrect(Form("%s/ZeeMC_PF_%s_2G/", directory.Data(), sqrts.Data()), i, rec_sig);
            rcStatW[i]->loadRooWorkspacesDiagData(Form("%s/ZeeData_PF_%s_2G/", directory.Data(), sqrts.Data()), i, rec_sig);
            rcStatW[i]->loadRooWorkspacesDiagMC(Form("%s/ZeeMC_PF_%s_2G/", directory.Data(), sqrts.Data()), i, rec_sig);
        }
    }
    if (doEta)
    {
        rcEta05Wp->loadRooWorkspacesMCtoCorrect(Form("%s/WepMC_PF_%s_Eta1/", directory.Data(), sqrts.Data()));
        rcEta05Wp->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Eta1/", directory.Data(), sqrts.Data()));
        rcEta05Wp->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Eta1/", directory.Data(), sqrts.Data()));

        rcEta051Wp->loadRooWorkspacesMCtoCorrect(Form("%s/WepMC_PF_%s_Eta2/", directory.Data(), sqrts.Data()));
        rcEta051Wp->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Eta2/", directory.Data(), sqrts.Data()));
        rcEta051Wp->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Eta2/", directory.Data(), sqrts.Data()));

        rcEta1Wp->loadRooWorkspacesMCtoCorrect(Form("%s/WepMC_PF_%s_Eta3/", directory.Data(), sqrts.Data()));
        rcEta1Wp->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Eta3/", directory.Data(), sqrts.Data()));
        rcEta1Wp->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Eta3/", directory.Data(), sqrts.Data()));

        rcEta05Wm->loadRooWorkspacesMCtoCorrect(Form("%s/WemMC_PF_%s_Eta1/", directory.Data(), sqrts.Data()));
        rcEta05Wm->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Eta1/", directory.Data(), sqrts.Data()));
        rcEta05Wm->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Eta1/", directory.Data(), sqrts.Data()));

        rcEta051Wm->loadRooWorkspacesMCtoCorrect(Form("%s/WemMC_PF_%s_Eta2/", directory.Data(), sqrts.Data()));
        rcEta051Wm->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Eta2/", directory.Data(), sqrts.Data()));
        rcEta051Wm->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Eta2/", directory.Data(), sqrts.Data()));

        rcEta1Wm->loadRooWorkspacesMCtoCorrect(Form("%s/WemMC_PF_%s_Eta3/", directory.Data(), sqrts.Data()));
        rcEta1Wm->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Eta3/", directory.Data(), sqrts.Data()));
        rcEta1Wm->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Eta3/", directory.Data(), sqrts.Data()));
    }
    if (doKeys)
    {
        rcKeysWp->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WepMC_PF_%s_Keys/", directory.Data(), sqrts.Data()));
        rcKeysWp->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Keys/", directory.Data(), sqrts.Data()));
        rcKeysWp->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Keys/", directory.Data(), sqrts.Data()));

        rcKeysWm->loadRooWorkspacesMCtoCorrectKeys(Form("%s/WemMC_PF_%s_Keys/", directory.Data(), sqrts.Data()));
        rcKeysWm->loadRooWorkspacesData(Form("%s/ZeeData_PF_%s_Keys/", directory.Data(), sqrts.Data()));
        rcKeysWm->loadRooWorkspacesMC(Form("%s/ZeeMC_PF_%s_Keys/", directory.Data(), sqrts.Data()));
    }

    // load the MET XY correction file
    // can the XY correction from Z's really can be applied to W's ?
    TString channel = "ee";
    METXYCorrector metXYCorr("XYCorrector", (envStr + "/MitEwk13TeV/Recoil/data/met_xy_" + sqrts + "_" + channel + ".root").Data());

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    // Create output directory
    gSystem->mkdir(outputDir, kTRUE);
    // CPlot::sOutDir = outputDir;

    // ----------------------------------------------------------------------------------------------------------------------------------

    TFile *infile = 0;
    TTree *intree = 0;

    // Read input file and get the TTrees
    cout << "Processing " << fileName.Data() << "..." << endl;
    infile = TFile::Open((inputDir + TString("/") + fileName).Data());
    assert(infile);
    intree = (TTree *)infile->Get("Events");
    assert(intree);

    TH1D *hGenWeights = (TH1D *)infile->Get("hGenWeights");
    TH1D *hLHEWeightSum = (TH1D *)infile->Get("hLHEWeightSum");

    // Variables to get some of the branches out of the tree
    Float_t scale1fb, scale1fbUp, scale1fbDown, prefireWeight, prefireUp, prefireDown;
    Float_t prefireEcal, prefireEcalUp, prefireEcalDown, prefireMuon, prefireMuonUp, prefireMuonDown;
    Float_t met, metPhi, sumEt, mt, u1, u2;
    Int_t q;
    TLorentzVector *lep = 0, *genV = 0, *genLep = 0;
    TLorentzVector *lep_raw = 0;
    Float_t pfCombIso;
    Float_t lepError; //, lep2error;

    intree->SetBranchAddress("prefireWeight", &prefireWeight);     // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireUp", &prefireUp);             // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireDown", &prefireDown);         // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcal", &prefireEcal);         // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcalUp", &prefireEcalUp);     // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcalDown", &prefireEcalDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuon", &prefireMuon);         // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuonUp", &prefireMuonUp);     // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuonDown", &prefireMuonDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fb", &scale1fb);               // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbUp", &scale1fbUp);           // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fbDown", &scale1fbDown);       // event weight per 1/fb (MC)
    intree->SetBranchAddress(met_name.c_str(), &met);              // MET
    intree->SetBranchAddress(metPhi_name.c_str(), &metPhi);        // phi(MET)
    // intree->SetBranchAddress(u1_name.c_str(), &u1); // parallel component of recoil
    // intree->SetBranchAddress(u2_name.c_str(), &u2); // perpendicular component of recoil
    intree->SetBranchAddress("q", &q);                 // lepton charge
    intree->SetBranchAddress("lep", &lep);             // lepton 4-vector
    intree->SetBranchAddress("lep_raw", &lep_raw);     // lepton 4-vector (raw)
    intree->SetBranchAddress("genLep", &genLep);       // lepton 4-vector
    intree->SetBranchAddress("genV", &genV);           // lepton 4-vector
    intree->SetBranchAddress("pfCombIso", &pfCombIso); // lepton 4-vector
    intree->SetBranchAddress("lepError", &lepError);   // sc2 4-vector

    Long64_t nevents = intree->GetEntries();
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

    //
    // Set up output file
    //
    TString outfilename = outputDir + TString("/") + fileName;
    if (NSEC != 1)
    {
        // file broken into different sections;
        outfilename.Remove(outfilename.Length() - 5, 5);
        outfilename += TString("_") + Form("%dSections", NSEC) + TString("_") + Form("%d", ITH) + TString(".root");
    }
    TFile *outFile = TFile::Open(outfilename, "RECREATE");
    TH1::AddDirectory(kFALSE);
    // Actually just need to clone the tree:
    TTree *outTree = intree->CloneTree(0);
    // Variables for the new branches:
    Double_t effSFweight = 1, relIso = 0;
    Double_t evtWeightSysFSR = 1, evtWeightSysMC = 1, evtWeightSysBkg = 1;
    Double_t mtCorr = 0;
    vector<Double_t> metVars, metVarsPhi;
    vector<Double_t> evtWeight;

    for (int i = 0; i < nMET; i++)
    {
        metVars.push_back(0);
        metVarsPhi.push_back(0);
    }
    for (int i = 0; i < nWeight; i++)
        evtWeight.push_back(0);

    outFile->cd();
    outTree->Branch("relIso", &relIso, "relIso/d", 99);                 // scaled isolation variable that needs calculation
    outTree->Branch("mtCorr", &mtCorr, "mtCorr/d", 99);                 // corrected MET with keys corrections
    outTree->Branch("evtWeight", "vector<Double_t>", &evtWeight, 99);   // event weight vector
    outTree->Branch("effSFweight", &effSFweight, "effSFweight/d", 99);  // scale factors weight
    outTree->Branch("lep_raw", "TLorentzVector", &lep_raw, 99);         // uncorrected lepton vector
    outTree->Branch("metVars", "vector<Double_t>", &metVars, 99);       // uncorrected lepton vector
    outTree->Branch("metVarsPhi", "vector<Double_t>", &metVarsPhi, 99); // uncorrected lepton vector

    // Double_t mt=-999;

    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << ", among these processing " << IEND - IBEGIN << std::endl;
    // for (UInt_t ientry = 0; ientry < intree->GetEntries(); ientry++) {
    for (Long64_t ientry = IBEGIN; ientry < IEND; ientry++)
    {
        intree->GetEntry(ientry);
        if (ientry % 10000 == 0)
            cout << "Event " << ientry << ". " << (double)ientry / (double)intree->GetEntries() * 100 << " % done with this file." << endl;

        // vector containing raw lepton info for correcting MET
        TVector2 vLepRaw((lep_raw->Pt()) * cos(lep_raw->Phi()), (lep_raw->Pt()) * sin(lep_raw->Phi()));

        double pU1 = 0; //--
        double pU2 = 0; //--

        // data/MC scale factor corrections
        Double_t effdata, effmc;
        Double_t effdatah, effmch, effdatal, effmcl;
        Double_t effdataFSR, effdataMC, effdataBkg;
        Double_t edTag, emTag;
        Double_t corr = 1, corrdu = 1, corrdd = 1, corrmu = 1, corrmd = 1;
        Double_t corrFSR = 1, corrMC = 1, corrBkg = 1, corrTag = 1;

        if (fabs(lep->Pt()) < PT_CUT)
            continue;

        if (fabs(lep->Eta()) > ETA_CUT)
            continue;

        if (fabs(lep->Eta()) >= ECAL_GAP_LOW && fabs(lep->Eta()) <= ECAL_GAP_HIGH)
            continue;

        // org is the original MET in the ntuple
        // without any correction
        metVars[org] = met;
        metVarsPhi[org] = metPhi;

        if (doMETXYCorrection)
        {
            // apply MET XY correction before any MET corrections
            // std::cout << "before correction met = " << met << " metPhi = " << metPhi << std::endl;
            metXYCorr.CorrectMETXY(met, metPhi, isData);
            // std::cout << "after correction met = " << met << " metPhi = " << metPhi << std::endl;
        }

        TVector2 vLepCor((lep->Pt()) * cos(lep->Phi()), (lep->Pt()) * sin(lep->Phi()));
        Double_t lepPt = vLepCor.Mod();
        TVector2 vMetCorr((met)*cos(metPhi), (met)*sin(metPhi));
        Double_t corrMet = (vMetCorr + vLepRaw - vLepCor).Mod();
        Double_t corrMetPhi = (vMetCorr + vLepRaw - vLepCor).Phi();

        if (isData)
        {
            mt = sqrt(2.0 * (lep->Pt()) * (corrMet) * (1.0 - cos(toolbox::deltaPhi(lep->Phi(), corrMetPhi))));

            metVars[no] = corrMet;
            metVars[ru] = corrMet;
            metVars[rd] = corrMet;
            metVarsPhi[no] = corrMetPhi;
            metVarsPhi[ru] = corrMetPhi;
            metVarsPhi[rd] = corrMetPhi;
            mtCorr = mt;
        }
        else
        {
            TLorentzVector vEle;
            vEle.SetPtEtaPhiM(lep->Pt(), lep->Eta(), lep->Phi(), ELE_MASS);

            corr = effs.fullCorrections(&vEle, q);
            vector<double> uncs_gsf = effs.getUncSel(&vEle, q);

            corrFSR *= uncs_gsf[0] * effs.computeHLTSF(&vEle, q); // alternate fsr model
            corrMC *= uncs_gsf[1] * effs.computeHLTSF(&vEle, q);  // alternate mc gen model
            corrBkg *= uncs_gsf[2] * effs.computeHLTSF(&vEle, q); // alternate bkg model
            corrTag *= uncs_gsf[3] * effs.computeHLTSF(&vEle, q); // alternate bkg model

            double var = 0.;
            var += effs.statUncHLT(&vEle, q, hErr, hErr, 1.0);
            var += effs.statUncSel(&vEle, q, hErr, hErr, 1.0, true);

            evtWeight[main] = corr * scale1fb * prefireWeight;
            evtWeight[fsr] = corrFSR * scale1fb * prefireWeight;
            evtWeight[mc] = corrMC * scale1fb * prefireWeight;
            evtWeight[bkg] = corrBkg * scale1fb * prefireWeight;
            evtWeight[tagpt] = corrTag * scale1fb * prefireWeight;
            evtWeight[effstat] = (corr + sqrt(fabs(var))) * scale1fb * prefireWeight; // not used
            evtWeight[pfireu] = corr * scale1fb * prefireUp;
            evtWeight[pfired] = corr * scale1fb * prefireDown;
            evtWeight[pfireecalu] = (prefireEcal > 0) ? (evtWeight[main] * prefireEcalUp / prefireEcal) : 0.;
            evtWeight[pfireecald] = (prefireEcal > 0) ? (evtWeight[main] * prefireEcalDown / prefireEcal) : 0.;
            evtWeight[pfiremuu] = (prefireMuon > 0) ? (evtWeight[main] * prefireMuonUp / prefireMuon) : 0.;
            evtWeight[pfiremud] = (prefireMuon > 0) ? (evtWeight[main] * prefireMuonDown / prefireMuon) : 0.;
            evtWeight[effstat_lepPos] = (q > 0) ? evtWeight[effstat] : evtWeight[main];
            evtWeight[effstat_lepNeg] = (q < 0) ? evtWeight[effstat] : evtWeight[main];

            TLorentzVector lepD, lepU;
            lepU.SetPtEtaPhiM(lep->Pt() * (1 + lepError), lep->Eta(), lep->Phi(), ELE_MASS);
            lepD.SetPtEtaPhiM(lep->Pt() * (1 - lepError), lep->Eta(), lep->Phi(), ELE_MASS);
            TVector2 vLepU(lepU.Pt() * cos(lepU.Phi()), lepU.Pt() * sin(lepU.Phi()));
            TVector2 vLepD(lepD.Pt() * cos(lepD.Phi()), lepD.Pt() * sin(lepD.Phi()));

            Double_t corrMetU = (vMetCorr + vLepRaw - vLepU).Mod();
            Double_t corrMetPhiU = (vMetCorr + vLepRaw - vLepU).Phi();
            Double_t corrMetD = (vMetCorr + vLepRaw - vLepD).Mod();
            Double_t corrMetPhiD = (vMetCorr + vLepRaw - vLepD).Phi();

            metVars[no] = corrMet;
            metVarsPhi[no] = corrMetPhi;
            metVars[keys] = corrMet;
            metVarsPhi[keys] = corrMetPhi;
            metVars[eta] = corrMet;
            metVarsPhi[eta] = corrMetPhi;
            metVars[cent] = corrMet;
            metVarsPhi[cent] = corrMetPhi;
            metVars[ru] = corrMetU;
            metVarsPhi[ru] = corrMetPhiU;
            metVars[rd] = corrMetD;
            metVarsPhi[rd] = corrMetPhiD;
            for (int i = 0; i < nNV; i++)
            {
                int ofs = i + ns;
                metVars[ofs] = corrMet;
                metVarsPhi[ofs] = corrMetPhi;
            }

            if (isRecoil)
            {
                double genVPt = genV->Pt();
                double genVPhi = genV->Phi();
                double genVy = genV->Rapidity();
                if (q > 0)
                {
                    if (doKeys)
                    {
                        rcKeysWp->CorrectInvCdf(metVars[keys], metVarsPhi[keys], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kTRUE, kFALSE);
                    }
                    if (doEta)
                    {
                        if (fabs(genVy) < 0.5)
                            rcEta05Wp->CorrectInvCdf(metVars[eta], metVarsPhi[eta], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        else if (fabs(genVy) >= 0.5 && fabs(genVy) < 1.0)
                            rcEta051Wp->CorrectInvCdf(metVars[eta], metVarsPhi[eta], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        else
                            rcEta1Wp->CorrectInvCdf(metVars[eta], metVarsPhi[eta], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                    }
                    if (doInclusive)
                    {
                        rcMainWp->CorrectInvCdf(metVars[cent], metVarsPhi[cent], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        rcMainWp->CorrectInvCdf(metVars[ru], metVarsPhi[ru], genVPt, genVPhi, lepU.Pt(), lepU.Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        rcMainWp->CorrectInvCdf(metVars[rd], metVarsPhi[rd], genVPt, genVPhi, lepD.Pt(), lepD.Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                    }
                    if (doStat)
                    {
                        for (int i = 0; i < nNV; i++)
                        {
                            int ofs = i + ns;
                            rcStatW[i]->CorrectInvCdf(metVars[ofs], metVarsPhi[ofs], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kTRUE);
                        }
                    }
                }
                else
                {
                    if (doKeys)
                    {
                        rcKeysWm->CorrectInvCdf(metVars[keys], metVarsPhi[keys], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kTRUE, kFALSE);
                    }
                    if (doEta)
                    {
                        // metVars[eta]=corrMet; metVarsPhi[eta]=corrMetPhi;
                        if (fabs(genVy) < 0.5)
                            rcEta05Wm->CorrectInvCdf(metVars[eta], metVarsPhi[eta], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        else if (fabs(genVy) >= 0.5 && fabs(genVy) < 1.0)
                            rcEta051Wm->CorrectInvCdf(metVars[eta], metVarsPhi[eta], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        else
                            rcEta1Wm->CorrectInvCdf(metVars[eta], metVarsPhi[eta], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                    }
                    if (doInclusive)
                    {
                        rcMainWm->CorrectInvCdf(metVars[cent], metVarsPhi[cent], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        rcMainWm->CorrectInvCdf(metVars[ru], metVarsPhi[ru], genVPt, genVPhi, lepU.Pt(), lepU.Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        rcMainWm->CorrectInvCdf(metVars[rd], metVarsPhi[rd], genVPt, genVPhi, lepD.Pt(), lepD.Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                        ;
                    }
                    if (doStat)
                    {
                        for (int i = 0; i < nNV; i++)
                        {
                            int ofs = i + ns;
                            rcStatW[i]->CorrectInvCdf(metVars[ofs], metVarsPhi[ofs], genVPt, genVPhi, lep->Pt(), lep->Phi(), pU1, pU2, 0, 0, 0, kFALSE, kTRUE);
                        }
                    }
                }
            }
            mtCorr = sqrt(2.0 * (lep->Pt()) * (metVars[cent]) * (1.0 - cos(toolbox::deltaPhi(lep->Phi(), metVarsPhi[cent]))));
        }

        // std::cout << "done w recoil" << std::endl;
        // set electron relIso to be these quantities but scaled to a constant 0.15 cutoff
        if (fabs(lep->Eta()) < ECAL_GAP_LOW)
            relIso = (pfCombIso - 0.506) / lep->Pt() + (0.15 - 0.0478);
        if (fabs(lep->Eta()) > ECAL_GAP_LOW)
            relIso = (pfCombIso - 0.963) / lep->Pt() + (0.15 - 0.0658);

        effSFweight = corr;
        outTree->Fill(); // add new info per event to the new tree
    }                    // end of loop over events
    delete rcMainWp;
    delete rcMainWm;
    delete rcKeysWp;
    delete rcKeysWm;
    for (int i = 0; i < nNV; i++)
        delete rcStatW[i];
    delete rcEta05Wp;
    delete rcEta051Wp;
    delete rcEta1Wp;
    delete rcEta05Wm;
    delete rcEta051Wm;
    delete rcEta1Wm;

    outFile->cd();
    hGenWeights->Write();
    if (hLHEWeightSum)
    {
        // only signal region contains hLHEWeightSum
        hLHEWeightSum->Write();
    }
    outFile->Write();
    std::cout << "wrote outfile" << std::endl;

    delete intree;
    delete infile;
    gBenchmark->Show("fitWe");
} // end of function
