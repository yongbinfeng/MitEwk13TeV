//================================================================================================
//
// Perform fit to extract Z->ee signal and efficiency simultaneously
//
//  * outputs plots and fit results summary
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
#include "MitEwk13TeV/EleScale/EnergyScaleCorrection.h"
#include "MitEwk13TeV/Utils/RecoilCorrector.hh"
#include "MitEwk13TeV/Utils/METXYCorrector.hh"

#include "MitEwk13TeV/Utils/AppEffSF.cc"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void ZeeNtupleMod(
    const TString outputDir, // output directory
    const TString inputDir, // input directory
    const TString sqrts,
    const TString fileName, // both the input and output final file name i.e. data_select.root
    const Int_t NSEC = 1,
    const Int_t ITH = 0
    )
{
    std::cout << "---------------- STARTING ZEE ------------------------" << std::endl;
    gBenchmark->Start("plotZee");

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
    const Double_t ELE_MASS = 0.000511;
    //
    // input ntuple file names
    //
    enum { eData,
        eZee,
        eEWK,
        eTop,
        eDib,
        eWx,
        eZxx }; // data type enum

    int filetype = -1;
    if (fileName.Contains("data_select")) {
        filetype = eData;
    } else if (fileName.Contains("zee_select")) {
        filetype = eZee;
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

    Bool_t isRecoil = (fileName.Contains("zee_select")  || fileName.Contains("wx_select")  || fileName.Contains("wx0_select")  || fileName.Contains("wx1_select")  || fileName.Contains("wx2_select")  || fileName.Contains("zxx_select"));
    std::cout << "do Recoil " << isRecoil << std::endl;

    const TString envStr = (TString)gSystem->Getenv("CMSSW_BASE") + "/src/";

    // apply the muon channel information for electrons
    // the assumption is electron channel is similar to muon channel
    const TString directory((envStr + "Corrections/Recoil").Data());
    RecoilCorrector* rcMainZ = new RecoilCorrector("", "");
    rcMainZ->loadRooWorkspacesMCtoCorrect(Form("%s/ZmmMC_PF_%s_2G/", directory.Data(), sqrts.Data()));
    rcMainZ->loadRooWorkspacesData(Form("%s/ZmmData_PF_%s_2G_bkg_fixRoch/", directory.Data(), sqrts.Data()));
    rcMainZ->loadRooWorkspacesMC(Form("%s/ZmmMC_PF_%s_2G/", directory.Data(), sqrts.Data()));

    //
    // warning: this might need to be updated
    //
    //METXYCorrector* metcorXY = new METXYCorrector("", "");
    //metcorXY->loadXYCorrection("/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_10_6_0/src/PostCorrNTuple/root/output_metxy.root");

    //
    // Fit options
    //
    const Double_t MASS_LOW = 60;
    const Double_t MASS_HIGH = 120;
    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4; //4;

    const Double_t ECAL_GAP_LOW  = 1.4442;
    const Double_t ECAL_GAP_HIGH = 1.566;
    //const Double_t ECAL_GAP_LOW = 10.;
    //const Double_t ECAL_GAP_HIGH = 10.;

    // efficiency files

    TString baseDir = envStr + "Corrections/Efficiency/" + sqrts + "/results/Zee/";
    AppEffSF effs(baseDir);
    effs.loadHLT("EleHLTEff_aMCxPythia", "Positive", "Negative");
    effs.loadSel("EleGSFSelEff_aMCxPythia", "Combined", "Combined");
    
    //
    // Warning: this would need to be updated for 5TeV
    //
    TString sysDir = envStr + "Corrections/Efficiency/" + sqrts + "/results/Systematics/";
    TString SysFileGSFSel = sysDir + "SysUnc_EleGSFSelEff.root";
    effs.loadUncSel(SysFileGSFSel);
    TH2D* hErr = new TH2D("hErr", "", 10, 0, 10, 20, 0, 20);

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    // event category enumeration
    enum { eEleEle2HLT = 1,
        eEleEle1HLT1L1,
        eEleEle1HLT,
        eEleEleNoSel,
        eEleSC };

    gSystem->mkdir(outputDir, kTRUE);

    //
    // Declare variables to read in ntuple
    //
    UInt_t runNum, lumiSec, evtNum;
    UInt_t matchGen;
    UInt_t category;
    UInt_t npv, npu;
    Float_t scale1fb;
    Float_t prefireWeight, prefireUp = 1, prefireDown = 1;
    Float_t prefireEcal, prefireEcalUp, prefireEcalDown, prefireMuon, prefireMuonUp, prefireMuonDown;
    Float_t met, metPhi, u1, u2;
    Float_t lep1error, lep2error;
    Int_t q1, q2;
    TLorentzVector *lep1 = 0, *lep2 = 0;
    TLorentzVector *lep1_raw = 0, *lep2_raw = 0;
    TLorentzVector *dilep = 0, *dilepSC = 0;
    TLorentzVector *genV = 0;
    TLorentzVector *sc1 = 0, *sc2 = 0;

    Float_t r91 = 0;
    Float_t r92 = 0;
    Float_t random = 0;

    Double_t nDib = 0, nWx = 0, nZxx = 0;
    Double_t nDibUnc = 0, nWxUnc = 0, nZxxUnc = 0;

    TFile* infile = 0;
    TTree* intree = 0;

    // Read input file and get the TTrees
    cout << "Processing " << fileName.Data() << "..." << endl;
    infile = TFile::Open((inputDir + TString("/") + fileName).Data());
    assert(infile);
    intree = (TTree*)infile->Get("Events");
    assert(intree);

    intree->SetBranchAddress("runNum", &runNum); // event run number
    intree->SetBranchAddress("lumiSec", &lumiSec); // event lumi section
    intree->SetBranchAddress("evtNum", &evtNum); // event number
    intree->SetBranchAddress("category", &category); // dilepton category
    intree->SetBranchAddress("npv", &npv); // number of primary vertices
    intree->SetBranchAddress("npu", &npu); // number of in-time PU events (MC)
    intree->SetBranchAddress("prefireWeight", &prefireWeight); // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("prefireUp", &prefireUp); // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("prefireDown", &prefireDown); // prefire weight for 2017 conditions (MC)
    intree->SetBranchAddress("prefireEcal", &prefireEcal); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcalUp", &prefireEcalUp); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireEcalDown", &prefireEcalDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuon", &prefireMuon); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuonUp", &prefireMuonUp); // event weight per 1/fb (MC)
    intree->SetBranchAddress("prefireMuonDown", &prefireMuonDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("scale1fb", &scale1fb); // event weight per 1/fb (MC)
    //intree->SetBranchAddress("scale1fbUp", &scale1fbUp); // event weight per 1/fb (MC)
    //intree->SetBranchAddress("scale1fbDown", &scale1fbDown); // event weight per 1/fb (MC)
    intree->SetBranchAddress("met", &met); // MET
    intree->SetBranchAddress("metPhi", &metPhi); // phi(MET)
    intree->SetBranchAddress("u1", &u1); // u1
    intree->SetBranchAddress("u2", &u2); // u2
    //intree->SetBranchAddress("genVMass", &genVMass); // event weight per 1/fb (MC)
    intree->SetBranchAddress("genV", &genV); // lepton 4-vector
    intree->SetBranchAddress("q1", &q1); // charge of tag lepton
    intree->SetBranchAddress("q2", &q2); // charge of probe lepton
    intree->SetBranchAddress("lep1_raw", &lep1); // tag lepton 4-vector
    intree->SetBranchAddress("lep2_raw", &lep2); // probe lepton 4-vector
    // intree->SetBranchAddress("lep1_raw",       &lep1_raw);        // tag lepton 4-vector
    // intree->SetBranchAddress("lep2_raw",       &lep2_raw);        // probe lepton 4-vector
    intree->SetBranchAddress("sc1", &sc1); // sc1 4-vector
    intree->SetBranchAddress("sc2", &sc2); // sc2 4-vector
    intree->SetBranchAddress("dilep", &dilep); // sc2 4-vector
    // intree->SetBranchAddress("dilepSC",   &dilepSC);      // sc2 4-vector
    intree->SetBranchAddress("r91", &r91); // sc2 4-vector
    intree->SetBranchAddress("r92", &r92); // sc2 4-vector
    // intree->SetBranchAddress("random",       &random);        // sc2 4-vector
    intree->SetBranchAddress("lep1error", &lep1error); // sc2 4-vector
    intree->SetBranchAddress("lep2error", &lep2error); // sc2 4-vector

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

    //
    // loop over events
    //
    std::cout << "Number of Events = " << intree->GetEntries() << ", among these processing " << IEND - IBEGIN << std::endl;
    //for (UInt_t ientry = 0; ientry < intree->GetEntries(); ientry++) {
    for (Long64_t ientry = IBEGIN; ientry < IEND; ientry++) {
        // for(UInt_t ientry=0; ientry<1000; ientry++) {
        if (ientry % 100000 == 0)
            cout << "Processing event " << ientry << ". " << (double)ientry / (double)intree->GetEntries() * 100 << " percent done with this file." << endl;
        intree->GetEntry(ientry);

        Double_t corr = 1;
        Double_t corrFSR = 1, corrMC = 1, corrBkg = 1, corrTag = 1;

        // std::cout << "r92 " << r92 << std::endl;
        if (fabs(lep1->Eta()) > ETA_CUT)
            continue;
        if (fabs(lep2->Eta()) > ETA_CUT)
            continue;

        if (q1 * q2 > 0)
            continue;
        if (lep1->Pt() < PT_CUT)
            continue;
        if (lep2->Pt() < PT_CUT)
            continue;

        if (isnan(prefireUp) || isnan(prefireDown)) {
            prefireUp = prefireWeight;
            prefireDown = prefireWeight;
        }

        if (fabs(lep1->Eta()) >= ECAL_GAP_LOW && fabs(lep1->Eta()) <= ECAL_GAP_HIGH)
            continue;
        if (fabs(lep2->Eta()) >= ECAL_GAP_LOW && fabs(lep2->Eta()) <= ECAL_GAP_HIGH)
            continue;

        if (!(category == 1) && !(category == 2) && !(category == 3))
            continue;
        if (filetype == eData) {

            TLorentzVector el1;
            TLorentzVector el2;
            el1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), ELE_MASS);
            el2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), ELE_MASS);

            Double_t lp1 = el1.Pt();
            Double_t lp2 = el2.Pt();
            Double_t lq1 = q1;
            Double_t lq2 = q2;

            TLorentzVector l1, l2;
            if (lp1 > lp2) {
                l1.SetPtEtaPhiM(lp1, lep1->Eta(), lep1->Phi(), ELE_MASS);
                l2.SetPtEtaPhiM(lp2, lep2->Eta(), lep2->Phi(), ELE_MASS);
            } else {
                l1.SetPtEtaPhiM(lp2, lep2->Eta(), lep2->Phi(), ELE_MASS);
                l2.SetPtEtaPhiM(lp1, lep1->Eta(), lep1->Phi(), ELE_MASS);
            }

            //TLorentzVector l1U, l2U;
            //l1U.SetPtEtaPhiM(lp1 * (1 + lep1error), lep1->Eta(), lep1->Phi(), ELE_MASS);
            //l2U.SetPtEtaPhiM(lp2 * (1 + lep2error), lep2->Eta(), lep2->Phi(), ELE_MASS);

            //TLorentzVector l1D, l2D;
            //l1D.SetPtEtaPhiM(lp1 * (1 - lep1error), lep1->Eta(), lep1->Phi(), ELE_MASS);
            //l2D.SetPtEtaPhiM(lp2 * (1 - lep2error), lep2->Eta(), lep2->Phi(), ELE_MASS);

            //double massU = (l1U + l2U).M();
            //double massD = (l1D + l2D).M();

            mass = (l1 + l2).M();
            // mass=dilepSC->M();
            
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
            //metcorXY->CalcU1U2(metVars[cxy], metVarsPhi[cxy], (l1+l2).Pt(), (l1+l2).Phi(), (l1+l2).Pt(), (l1+l2).Phi(), pU1_postXY, pU2_postXY);
        } else {
            Double_t lp1 = lep1->Pt();
            Double_t lp2 = lep2->Pt();
            Double_t lq1 = q1;
            Double_t lq2 = q2;
            // double rand1;

            // // set the smearings here
            // for(int i = 0; i < 5; ++i){
            // double rand = gRandom->Gaus(0,1);
            // if(i==2) rand1=rand;
            // }
            // hGausRandHere->Fill(rand1);
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), ELE_MASS);
            l2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), ELE_MASS);

            double mll = (l1 + l2).M();
            //if (mll < MASS_LOW)
            //    continue;
            //if (mll > MASS_HIGH)
            //    continue;
            // if(lp1        < PT_CUT)    continue;
            // if(lp2        < PT_CUT)    continue;

            // if(genVMass>MASS_LOW && genVMass<MASS_HIGH) continue;

            corr = effs.fullCorrections(&l1, q1, &l2, q2);
            // cout << corr1 << " " << corr << endl;

            vector<double> uncs_gsf = effs.getUncSel(&l1, q1, &l2, q2);

            corrFSR *= uncs_gsf[0] * effs.computeHLTSF(&l1, q1, &l2, q2); // alternate fsr model
            corrMC *= uncs_gsf[1] * effs.computeHLTSF(&l1, q1, &l2, q2); // alternate mc gen model
            corrBkg *= uncs_gsf[2] * effs.computeHLTSF(&l1, q1, &l2, q2); // alternate bkg model
            corrTag *= uncs_gsf[3] * effs.computeHLTSF(&l1, q1, &l2, q2); // alternate tag-pt model

            double var = 0.;
            var += effs.statUncSel(&l1, q1, hErr, hErr, 1.0, true);
            var += effs.statUncSel(&l2, q2, hErr, hErr, 1.0, true);
            //var += effs.statUncHLTDilep(&l1, q1, &l2, q2);
            var += effs.statUncHLT(&l1, q1, hErr, hErr, 1.0);
            var += effs.statUncHLT(&l2, q2, hErr, hErr, 1.0);

            double var_lep1 = 0.;
            var_lep1 += effs.statUncSel(&l1, q1, hErr, hErr, 1.0, true);
            var_lep1 += effs.statUncHLT(&l1, q1, hErr, hErr, 1.0);

            double var_lep2 = 0.;
            var_lep2 += effs.statUncSel(&l2, q2, hErr, hErr, 1.0, true);
            var_lep2 += effs.statUncHLT(&l2, q2, hErr, hErr, 1.0);

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

            // corr=1;
            mass = (l1 + l2).M();

            // recoil corrections
            metVars[no] = met;
            metVarsPhi[no] = metPhi;
            metVars[cent] = met;
            metVarsPhi[cent] = metPhi;

            // to save the corrected U1 and U2
            pU1 = u1;
            pU2 = u2;
            if (isRecoil) {
                rcMainZ->CorrectInvCdf(metVars[cent], metVarsPhi[cent], genV->Pt(), genV->Phi(), (l1+l2).Pt(), (l1+l2).Phi(), pU1, pU2, 0, 0, 0, kFALSE, kFALSE);
                //std::cout << "PU1 " << pU1 << " pU2 " << pU2 << std::endl;
            }

            // XY correction on MC
            metVars[cxy] =  metVars[cent];
            metVarsPhi[cxy] = metVarsPhi[cent];
            //metcorXY->CorrectMETXY(metVars[cxy], metVarsPhi[cxy], npv, 0);
            //metcorXY->CalcU1U2(metVars[cxy], metVarsPhi[cxy], (l1+l2).Pt(), (l1+l2).Phi(), genV->Pt(), genV->Phi(), pU1_postXY, pU2_postXY);
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

    gBenchmark->Show("plotZee");
}

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
