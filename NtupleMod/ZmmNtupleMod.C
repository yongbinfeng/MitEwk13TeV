//================================================================================================
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TLorentzVector.h" // 4-vector class
#include <TBenchmark.h> // class to track macro running statistics
#include <TFile.h> // file handle class
#include <TGaxis.h>
#include <TH1D.h> // histogram class
#include <TROOT.h> // access to gROOT, entry point to ROOT system
#include <TRandom3.h>
#include <TStyle.h> // class to handle ROOT plotting styles
#include <TSystem.h> // interface to OS
#include <TTree.h> // class to access ntuples
#include <fstream> // functions for file I/O
#include <iomanip> // functions to format standard I/O
#include <iostream> // standard I/O
#include <sstream> // class for parsing strings
#include <string> // C++ string class
#include <vector> // STL vector class

#include "MitEwk13TeV/Utils/LeptonCorr.hh" // Scale and resolution corrections
#include "MitEwk13TeV/Utils/MyTools.hh" // various helper functions
#include "MitEwk13TeV/Utils/RecoilCorrector.hh"
#include "MitEwk13TeV/Utils/METXYCorrector.hh"

#include "MitEwk13TeV/RochesterCorr/RoccoR.cc"
#include "MitEwk13TeV/Utils/AppEffSF.cc"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// make data-fit difference plots
// TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

void orderByLepPt(TLorentzVector& mu1, TLorentzVector& mu2, Int_t& q1, Int_t& q2, double mu_MASS);

//=== MAIN MACRO =================================================================================================

void ZmmNtupleMod(
    const TString outputDir, // output directory
    const TString inputDir, // input directory
    const TString sqrts,
    const TString fileName, // both the input and output final file name i.e. data_select.root
    const Int_t NSEC = 1,
    const Int_t ITH = 0
    )
{
    gBenchmark->Start("plotZmm");
    gStyle->SetTitleOffset(1.100, "Y");

    // Control the types of uncertainties
    enum { 
        no,
        cent,
        cxy
        };
    const string vMET[] = { "no", "cent", "cxy"};
    int nMET = sizeof(vMET) / sizeof(vMET[0]);

    enum { main,
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
        effstat_lepNeg};
    const string vWeight[] = { "eff", "mc", "fsr", "bkg", "tagpt", "effstat", "pfireu", "pfired", "pfireecalu", "pfireecald", "pfiremuu", "pfiremud", "effstat_lepPos", "effstat_lepNeg"};
    int nWeight = sizeof(vWeight) / sizeof(vWeight[0]);

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================
    const Double_t mu_MASS = 0.1057;
    //
    // input ntuple file names
    //
    enum { eData,
        eZmm,
        eEWK,
        eTop,
        eDib,
        eZxx,
        eWx }; // data type enum

    int filetype = -1;
    if (fileName.Contains("data_select")) {
        filetype = eData;
    } else if (fileName.Contains("zmm_select")) {
        filetype = eZmm;
    } else if (fileName.Contains("top_select") || fileName.Contains("top1_select") || fileName.Contains("top2_select") || fileName.Contains("top3_select")) {
        filetype = eTop;
    } else if (fileName.Contains("zz_select") || fileName.Contains("ww_select") || fileName.Contains("wz_select")) {
        filetype = eDib;
    } else if (fileName.Contains("zxx_select")) {
        filetype = eZxx;
    } else if (fileName.Contains("wx_select") || fileName.Contains("wx0_select") || fileName.Contains("wx1_select") || fileName.Contains("wx2_select")) {
        filetype = eWx;
    }
    std::cout << "data type " << filetype << std::endl;

    Bool_t isRecoil = (fileName.Contains("zmm_select") || fileName.Contains("wx_select")  || fileName.Contains("wx0_select")  || fileName.Contains("wx1_select")  || fileName.Contains("wx2_select")  || fileName.Contains("zxx_select")  || fileName.Contains("zmm_select")  || fileName.Contains("wx_select")  || fileName.Contains("wx0_select")  || fileName.Contains("wx1_select")  || fileName.Contains("wx2_select")  || fileName.Contains("zxx_select") );
    std::cout << "do Recoil " << isRecoil << std::endl;

    //
    // Fit options
    //
    const Double_t MASS_LOW = 60;
    const Double_t MASS_HIGH = 120;
    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;

    const bool doRoch = true;
    const TString envStr = (TString)gSystem->Getenv("CMSSW_BASE") + "/src/";

    // efficiency files
    TString baseDir = envStr + "Corrections/Efficiency/" + sqrts + "/results/Zmm/";
    AppEffSF effs(baseDir);
    effs.loadHLT("MuHLTEff_aMCxPythia", "Positive", "Negative");
    effs.loadSel("MuSITEff_aMCxPythia", "Combined", "Combined");
    effs.loadSta("MuStaEff_aMCxPythia", "Combined", "Combined");

    TString sysDir = envStr + "Corrections/Efficiency/" + sqrts + "/results/Systematics/";
    TString sysFileSIT = sysDir + "SysUnc_MuSITEff.root";
    TString sysFileSta = sysDir + "SysUnc_MuStaEff.root";
    effs.loadUncSel(sysFileSIT);
    effs.loadUncSta(sysFileSta);
    TH2D* hErr = new TH2D("hErr", "", 10, 0, 10, 20, 0, 20);
    
    const TString directory((envStr + "Corrections/Recoil").Data());

    RecoilCorrector* rcMainZ = new RecoilCorrector("", "");
    rcMainZ->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMC_PF_%s_2G/", directory.Data(), sqrts.Data()));
    rcMainZ->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G_bkg_fixRoch/", directory.Data(), sqrts.Data()));
    rcMainZ->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_2G/", directory.Data(), sqrts.Data()));

    //
    // Warning: this might need to be updated for 5TeV
    //
    //METXYCorrector* metcorXY = new METXYCorrector("", "");
    //metcorXY->loadXYCorrection("/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_10_6_0/src/PostCorrNTuple/root/output_metxy.root");


    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    // event category enumeration
    enum { eMuMu2HLT = 1,
        eMuMu1HLT1mu1,
        eMuMu1HLT,
        eMuMuNoSel,
        eMuSta,
        eMuTrk }; // event category enum

    // // histograms for full selection
    // double ZPtBins[]={0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,400,1000};
    // double Lep1PtBins[]={25,26.3,27.6,28.9,30.4,31.9,33.5,35.2,36.9,38.8,40.7,42.8,44.9,47.1,49.5,52.0,54.6,57.3,60.7,65.6,72.2,80.8,92.1,107,126,150,200,300};

    // Create output directory
    gSystem->mkdir(outputDir, kTRUE);

    //
    // Declare variables to read in ntuple
    //
    UInt_t runNum, lumiSec, evtNum;
    UInt_t matchGen;
    UInt_t category;
    Float_t scale1fb;
    Float_t prefireWeight, prefireUp, prefireDown;
    Float_t prefireEcal, prefireEcalUp, prefireEcalDown, prefireMuon, prefireMuonUp, prefireMuonDown;
    Float_t met, metPhi, u1, u2;
    Int_t q1, q2;
    UInt_t nTkLayers1, nTkLayers2;
    TLorentzVector *lep1 = 0, *lep2 = 0;
    TLorentzVector *genlep1 = 0, *genlep2 = 0;
    TLorentzVector *genV = 0;
    Float_t genMuonPt1, genMuonPt2;

    Double_t nDib = 0, nWx = 0, nZxx = 0;
    Double_t nDibUnc = 0, nWxUnc = 0, nZxxUnc = 0;
    UInt_t npv;

    // Loading the Rochster Corrections
    RoccoR rc((envStr + "/MitEwk13TeV/RochesterCorr/RoccoR2017.txt").Data());

    TFile* infile = 0;
    TTree* intree = 0;

    // Read input file and get the TTrees
    cout << "Processing " << fileName.Data() << "..." << endl;
    infile = TFile::Open((inputDir + TString("/") + fileName).Data());
    assert(infile);
    intree = (TTree*)infile->Get("Events");
    assert(intree);

    intree->SetBranchAddress("category", &category); // dilepton category
    intree->SetBranchAddress("prefireWeight", &prefireWeight); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireUp", &prefireUp); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireDown", &prefireDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcal", &prefireEcal); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcalUp", &prefireEcalUp); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcalDown", &prefireEcalDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuon", &prefireMuon); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuonUp", &prefireMuonUp); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuonDown", &prefireMuonDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fb", &scale1fb); // event weight per 1/fb (MC)
    intree->SetBranchAddress("met", &met); // MET
    intree->SetBranchAddress("metPhi", &metPhi); // phi(MET)
    intree->SetBranchAddress("u1", &u1); // u1
    intree->SetBranchAddress("u2", &u2); // u2
    intree->SetBranchAddress("q1", &q1); // charge of tag lepton
    intree->SetBranchAddress("q2", &q2); // charge of probe lepton
    intree->SetBranchAddress("lep1", &lep1); // tag lepton 4-vector
    intree->SetBranchAddress("lep2", &lep2); // probe lepton 4-vector
    intree->SetBranchAddress("genlep1", &genlep1); // tag lepton 4-vector
    intree->SetBranchAddress("genlep2", &genlep2); // probe lepton 4-vector
    intree->SetBranchAddress("genMuonPt1", &genMuonPt1); // probe lepton 4-vector
    intree->SetBranchAddress("genMuonPt2", &genMuonPt2); // probe lepton 4-vector
    intree->SetBranchAddress("genV", &genV); // lepton 4-vector
    intree->SetBranchAddress("nTkLayers1", &nTkLayers1);
    intree->SetBranchAddress("nTkLayers2", &nTkLayers2);
    intree->SetBranchAddress("npv", &npv);

    TH1D* hGenWeights;
    hGenWeights = (TH1D*)infile->Get("hGenWeights");

    Long64_t nevents = intree->GetEntries();
    Long64_t IBEGIN = 0;
    Long64_t IEND = nevents;

    if (NSEC!=1) {
        cout << "n sections " << NSEC << " ith part " << ITH << endl;
        Long64_t nevents_per_job = nevents / NSEC;
        Long64_t remainder = nevents % NSEC;
        IBEGIN = ITH * nevents_per_job;
        IEND = (ITH + 1) * nevents_per_job;
        if (ITH == NSEC - 1) {
            // for the last section, add the remaining into these;
            IEND = nevents;
        }
        cout << "start, end events: " << IBEGIN << " " << IEND << endl;
    }

    //
    // Set up output file
    //
    TString outfilename = outputDir + TString("/") + fileName;
    if (NSEC!=1) {
        // file broken into different sections;
        outfilename.Remove(outfilename.Length()-5, 5);
        outfilename += TString("_") + Form("%dSections", NSEC) + TString("_") + Form("%d", ITH) + TString(".root");
    }
    TFile* outFile = TFile::Open(outfilename, "RECREATE");
    TH1::AddDirectory(kFALSE);
    // Actually just need to clone the tree:
    TTree* outTree = intree->CloneTree(0);
    Double_t mass = 0;
    Double_t effSFweight = 1;
    Double_t pU1 = 0, pU2 = 0, pU1_postXY = 0, pU2_postXY = 0;
    Double_t lep1error, lep2error;
    vector<Double_t> evtWeight;
    vector<Double_t> metVars, metVarsPhi;

    for (int i = 0; i < nMET; i++) {
        metVars.push_back(0);
        metVarsPhi.push_back(0);
    }

    for (int i = 0; i < nWeight; i++)
        evtWeight.push_back(0);

    outFile->cd();
    outTree->Branch("mass", &mass, "mass/d", 99); // invariant mass of the two leptons 
    outTree->Branch("evtWeight", "vector<Double_t>", &evtWeight, 99); // event weight vector
    outTree->Branch("effSFweight", &effSFweight, "effSFweight/d", 99); // scale factors weight
    outTree->Branch("metVars", "vector<Double_t>", &metVars, 99); // metVars
    outTree->Branch("metVarsPhi", "vector<Double_t>", &metVarsPhi, 99); // metVars in phi
    outTree->Branch("pU1", &pU1, "pU1/d", 99); // u1 after recoil correction
    outTree->Branch("pU2", &pU2, "pU2/d", 99); // u2 after recoil correction
    outTree->Branch("pU1_postXY", &pU1_postXY, "pU1_postXY/d", 99); // u1 after recoil and metxy correction
    outTree->Branch("pU2_postXY", &pU2_postXY, "pU2_postXY/d", 99); // u2 after recoil and metxy correction
    outTree->Branch("lep1error", &lep1error, "lep1error/d", 99);
    outTree->Branch("lep2error", &lep2error, "lep2error/d", 99);

    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << ", among these processing " << IEND - IBEGIN << std::endl;
    for (Long64_t ientry = IBEGIN; ientry < IEND; ientry++) {
    //for (UInt_t ientry = 0; ientry < intree->GetEntries(); ientry++) {
        if (ientry % 100000 == 0)
            cout << "Processing event " << ientry << ". " << (double)ientry / (double)intree->GetEntries() * 100 << " percent done with this file." << endl;
        intree->GetEntry(ientry);

        Double_t corr = 1;
        Double_t corrFSR = 1, corrMC = 1, corrBkg = 1, corrTag = 1;

        lep1error = 0;
        lep2error = 0;

        // fill Z events passing selection
        if (!(category == eMuMu2HLT) && !(category == eMuMu1HLT) && !(category == eMuMu1HLT1mu1))
            continue;
        // cout << "pass trigger?" << endl;

        if (filetype == eData) {

            TLorentzVector mu1;
            TLorentzVector mu2;
            mu1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), mu_MASS);
            mu2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), mu_MASS);

            double dtSF1 = rc.kScaleDT(q1, mu1.Pt(), mu1.Eta(), mu1.Phi()); //, s=0, m=0);
            double dtSF2 = rc.kScaleDT(q2, mu2.Pt(), mu2.Eta(), mu2.Phi()); //s=0, m=0);
            if (doRoch) {
                mu1 *= dtSF1;
                mu2 *= dtSF2;
                (*lep1) *= dtSF1;
                (*lep2) *= dtSF2; 
                lep1error = rc.kScaleDTerror(q1, mu1.Pt(), mu1.Eta(), mu1.Phi());
                lep2error = rc.kScaleDTerror(q2, mu2.Pt(), mu2.Eta(), mu2.Phi());
            }

            if (mu1.Pt() < PT_CUT) 
                continue;
            if (mu2.Pt() < PT_CUT)
                continue;
            if (q1 * q2 > 0)
                continue;

            Double_t lp1 = mu1.Pt();
            Double_t lp2 = mu2.Pt();
            Int_t q1 = q1;
            Int_t q2 = q2;

            mass = (mu1 + mu2).M();

            // skip the impactos of lepton energy scale correction on MET as they are minor
            metVars[no] = met;
            metVarsPhi[no] = metPhi;

            // for data, no correction on recoil
            metVars[cent] = met;
            metVarsPhi[cent] = metPhi;
            pU1 = u1;
            pU2 = u2;

            // correct MET XY for data
            metVars[cxy] = met;
            metVarsPhi[cxy] = metPhi;
            //metcorXY->CorrectMETXY(metVars[cxy], metVarsPhi[cxy], npv, 1);
            //metcorXY->CalcU1U2(metVars[cxy], metVarsPhi[cxy], (mu1+mu2).Pt(), (mu1+mu2).Phi(), (mu1+mu2).Pt(), (mu1+mu2).Phi(), pU1_postXY, pU2_postXY);

        } else {

            TLorentzVector mu1;
            TLorentzVector mu2;
            mu1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), mu_MASS);
            mu2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), mu_MASS);
            TLorentzVector mu1u, mu1d, mu2u, mu2d;
            mu1u = mu1;
            mu1d = mu1;
            mu2u = mu2;
            mu2d = mu2;


            double mcSF1 = 1;
            double mcSF2 = 1;
            if (genMuonPt1 > 0) {
                mcSF1 = rc.kSpreadMC(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt1);
                lep1error = rc.kSpreadMCerror(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt1);
            } else {
                double rand = gRandom->Uniform(1);
                mcSF1 = rc.kSmearMC(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), nTkLayers1, rand);
                lep1error = rc.kSmearMCerror(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), nTkLayers1, rand);
            }
            if (genMuonPt2 > 0) {
                mcSF2 = rc.kSpreadMC(q2, mu2.Pt(), mu2.Eta(), mu2.Phi(), genMuonPt2);
                lep2error = rc.kSpreadMCerror(q2, mu2.Pt(), mu2.Eta(), mu2.Phi(), genMuonPt2);
            } else {
                double rand = gRandom->Uniform(1);
                mcSF2 = rc.kSmearMC(q2, mu2.Pt(), mu2.Eta(), mu2.Phi(), nTkLayers2, rand);
                lep2error = rc.kSmearMCerror(q2, mu2.Pt(), mu2.Eta(), mu2.Phi(), nTkLayers2, rand);
            }
            mu1 *= mcSF1;
            (*lep1) *= mcSF1;
            mu2 *= mcSF2;
            (*lep2) *= mcSF2;

            if (mu1.Pt() < PT_CUT)
                continue;
            if (mu2.Pt() < PT_CUT)
                continue;
            if (q1 * q2 > 0)
                continue;

            //double deltaMcSF1 = rc.kSpreadMCerror(q1, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt1);
            //double deltaMcSF2 = rc.kSpreadMCerror(q2, mu1.Pt(), mu1.Eta(), mu1.Phi(), genMuonPt2);
            //mu1u *= (1 + deltaMcSF1);
            //mu1d *= (1 - deltaMcSF1);
            //mu2u *= (1 + deltaMcSF2);
            //mu2d *= (1 - deltaMcSF2);

            //Double_t lp1 = mu1.Pt();
            //Double_t lp2 = mu2.Pt();

            //double massU = (mu1u + mu2u).M();
            //double massD = (mu1d + mu2d).M();

            double mll = (mu1 + mu2).M();
            corr = effs.fullCorrections(&mu1, q1, &mu2, q2);

            vector<double> uncs_sta = effs.getUncSta(&mu1, q1, &mu2, q2);
            vector<double> uncs_sit = effs.getUncSel(&mu1, q1, &mu2, q2);


            corrFSR *= uncs_sta[0] * uncs_sit[0] * effs.computeHLTSF(&mu1, q1, &mu2, q2); // alternate fsr model
            corrMC *= uncs_sta[1] * uncs_sit[1] * effs.computeHLTSF(&mu1, q1, &mu2, q2); // alternate mc gen model
            corrBkg *= uncs_sta[2] * uncs_sit[2] * effs.computeHLTSF(&mu1, q1, &mu2, q2); // alternate bkg model
            corrTag *= uncs_sta[3] * uncs_sit[3] * effs.computeHLTSF(&mu1, q1, &mu2, q2); // alternate bkg model


            double var = 0.;
            var += effs.statUncSta(&mu1, q1, hErr, hErr, 1.0, true);
            var += effs.statUncSta(&mu2, q2, hErr, hErr, 1.0, true);
            var += effs.statUncSel(&mu1, q1, hErr, hErr, 1.0, true);
            var += effs.statUncSel(&mu2, q2, hErr, hErr, 1.0, true);
            //var += effs.statUncHLTDilep(&mu1, q1, &mu2, q2);
            var += effs.statUncHLT(&mu1, q1, hErr, hErr, 1.0);
            var += effs.statUncHLT(&mu2, q2, hErr, hErr, 1.0);

            double var_lep1 = 0.;
            var_lep1 += effs.statUncSta(&mu1, q1, hErr, hErr, 1.0, true);
            var_lep1 += effs.statUncSel(&mu1, q1, hErr, hErr, 1.0, true);
            var_lep1 += effs.statUncHLT(&mu1, q1, hErr, hErr, 1.0);

            double var_lep2 = 0.;
            var_lep2 += effs.statUncSta(&mu2, q2, hErr, hErr, 1.0, true);
            var_lep2 += effs.statUncSel(&mu2, q2, hErr, hErr, 1.0, true);
            var_lep2 += effs.statUncHLT(&mu2, q2, hErr, hErr, 1.0);

            //std::cout << "stational stat unc " << effs.statUncSta(&mu1, q1, hErr, hErr, 1.0) << " and " << effs.statUncSta(&mu2, q2, hErr, hErr, 1.0) << std::endl;
            //std::cout << "selection stat unc " << effs.statUncSel(&mu1, q1, hErr, hErr, 1.0) << " and " << effs.statUncSel(&mu2, q2, hErr, hErr, 1.0) << std::endl;
            //std::cout << "hlt stat unc " << effs.statUncHLTDilep(&mu1, q1, &mu2, q2) << " lep1 " << effs.statUncHLT(&mu1, q1, hErr, hErr, 1.0) << " lep2 " << effs.statUncHLT(&mu2, q2, hErr, hErr, 1.0) << std::endl;

            evtWeight[main] = corr * scale1fb * prefireWeight;
            evtWeight[fsr] = corrFSR * scale1fb * prefireWeight;
            evtWeight[mc] = corrMC * scale1fb * prefireWeight;
            evtWeight[bkg] = corrBkg * scale1fb * prefireWeight;
            evtWeight[tagpt] = corrTag * scale1fb * prefireWeight;
            evtWeight[effstat] = (corr + sqrt(abs(var))) * scale1fb * prefireWeight;
            evtWeight[pfireu] = corr * scale1fb * prefireUp;
            evtWeight[pfired] = corr * scale1fb * prefireDown;
            evtWeight[pfireecalu] = (prefireEcal > 0) ? (evtWeight[main] * prefireEcalUp / prefireEcal) : 0.;
            evtWeight[pfireecald] = (prefireEcal > 0) ? (evtWeight[main] * prefireEcalDown / prefireEcal) : 0.;
            evtWeight[pfiremuu] = (prefireMuon > 0) ? (evtWeight[main] * prefireMuonUp / prefireMuon) : 0.;
            evtWeight[pfiremud] = (prefireMuon > 0) ? (evtWeight[main] * prefireMuonDown / prefireMuon) : 0.;

            double corr_lep1 = (corr + sqrt(abs(var_lep1))) * scale1fb * prefireWeight;
            double corr_lep2 = (corr + sqrt(abs(var_lep2))) * scale1fb * prefireWeight;
            evtWeight[effstat_lepPos] = (q1 > 0) ? corr_lep1 : corr_lep2;
            evtWeight[effstat_lepNeg] = (q1 < 0) ? corr_lep1 : corr_lep2; 

            mass = (mu1 + mu2).M();

            // recoil corrections
            metVars[no] = met;
            metVarsPhi[no] = metPhi;
            metVars[cent] = met;
            metVarsPhi[cent] = metPhi;

            // to save the corrected U1 and U2
            pU1 = u1;
            pU2 = u2;
            if (isRecoil) {
                rcMainZ->CorrectInvCdf(metVars[cent], metVarsPhi[cent], genV->Pt(), genV->Phi(), (mu1+mu2).Pt(), (mu1+mu2).Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                //std::cout << "PU1 " << pU1 << " pU2 " << pU2 << std::endl;
            }

            // XY correction on MC
            metVars[cxy] =  metVars[cent];
            metVarsPhi[cxy] = metVarsPhi[cent];
            //metcorXY->CorrectMETXY(metVars[cxy], metVarsPhi[cxy], npv, 0);
            //metcorXY->CalcU1U2(metVars[cxy], metVarsPhi[cxy], (mu1+mu2).Pt(), (mu1+mu2).Phi(), genV->Pt(), genV->Phi(), pU1_postXY, pU2_postXY);
        }
        effSFweight = corr;
        outTree->Fill(); // add new info per event to the new tree
    }

    outFile->cd();
    hGenWeights->Write();
    outFile->Write();
    std::cout << "wrote outfile" << std::endl;

    delete infile;
    infile = 0, intree = 0;

    gBenchmark->Show("plotZmm");
}

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

void orderByLepPt(TLorentzVector& mu1, TLorentzVector& mu2, Int_t& q1, Int_t& q2, double mu_MASS)
{
    if (mu2.Pt() > mu1.Pt()) {
        TLorentzVector tmp;
        tmp.SetPtEtaPhiM(mu1.Pt(), mu1.Eta(), mu1.Phi(), mu_MASS);
        mu1 = mu2;
        mu2 = tmp;

        q1 *= -1;
        q2 *= -1;
    }
    return;
}
