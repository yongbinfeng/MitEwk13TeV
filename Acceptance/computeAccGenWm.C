//================================================================================================
// Not used for 13 TeV measurement.
//
// Compute W->munu acceptance at generator level
//
//  * outputs results summary text file
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TBenchmark.h> // class to track macro running statistics
#include <TCanvas.h>
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

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "MitEwk13TeV/Utils/CSample.hh" // helper class to handle samples
#include "MitEwk13TeV/Utils/ConfParse.hh" // input conf file parser
#include "MitEwk13TeV/Utils/MyTools.hh"
#endif

//=== MAIN MACRO =================================================================================================

void computeAccGenWm(const TString conf, // input file
    const TString outputDir, // output directory
    const TString outputName, // output filename
    const bool doDressed = 0,
    const Int_t charge = 0, // 0 = inclusive, +1 = W+, -1 = W-
    const bool doborn = 0)
{
    gBenchmark->Start("computeAccGenWm");

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================

    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;
    const Double_t ETA_BARREL = 1.2;
    const Double_t ETA_ENDCAP = 1.2;

    const Int_t BOSON_ID = 24;
    const Int_t LEPTON_ID = 13;

    const Int_t NPDF = 100;
    const Int_t NQCD = 6;

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
    baconhep::TGenEventInfo* gen = new baconhep::TGenEventInfo();
    TClonesArray* genPartArr = new TClonesArray("baconhep::TGenParticle");

    TFile* infile = 0;
    TTree* eventTree = 0;

    Double_t nEvtsv = 0, nSelv = 0, nSelBv = 0, nSelEv = 0;
    Double_t nEntries = 0, nEvtsAfter1Lep = 0, nEvtsAfterMT = 0;

    Double_t accv = 0, accBv = 0, accEv = 0;
    Double_t accErrv = 0, accErrBv = 0, accErrEv = 0;

    Double_t nEvtsv_pT = 0, nSelv_pT = 0, accv_pT = 0;

    vector<Double_t> nEvtsv_QCD(NQCD, 0), nSelv_QCD(NQCD, 0);
    vector<Double_t> nEvtsv_PDF(NPDF, 0), nSelv_PDF(NPDF, 0);

    TString sqrts = "13TeV";
    if (conf.Contains("5"))
        sqrts = "5TeV";
    //TFile* rf = new TFile("/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction/Z_pT/zPt_Normal" + sqrts + ".root");
    //TH1D* hh_diff = (TH1D*)rf->Get("hZptRatio");

    CSample* samp = 0;

    //
    // loop over the processes
    //
    for (UInt_t isamp = 0; isamp < samplev.size(); ++isamp) {

        //
        // loop through files
        //
        samp = samplev[isamp];
        const UInt_t nfiles = samp->fnamev.size();

        for (UInt_t ifile = 0; ifile < nfiles; ifile++) {

            // Read input file and get the TTrees
            cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... ";
            cout.flush();

            infile = TFile::Open(samp->fnamev[ifile]);
            assert(infile);

            eventTree = (TTree*)infile->Get("Events");
            assert(eventTree);
            eventTree->SetBranchAddress("GenEvtInfo", &gen);
            TBranch* genBr = eventTree->GetBranch("GenEvtInfo");
            eventTree->SetBranchAddress("GenParticle", &genPartArr);
            TBranch* partBr = eventTree->GetBranch("GenParticle");

            //
            // loop over events
            //
            double frac = 0.05; // fraction of events to be used for calculation
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
                genPartArr->Clear();
                partBr->GetEntry(ientry);

                Int_t lepq1 = -99;
                Int_t lepq2 = -99;
                if (fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                    continue;
                if (charge == -1 && toolbox::flavor(genPartArr, BOSON_ID) != LEPTON_ID)
                    continue;
                if (charge == 1 && toolbox::flavor(genPartArr, BOSON_ID) != -LEPTON_ID)
                    continue;
                if (charge == 0 && fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                    continue;

                TLorentzVector* vec = new TLorentzVector(0, 0, 0, 0);
                TLorentzVector* lep1 = new TLorentzVector(0, 0, 0, 0);
                TLorentzVector* lep2 = new TLorentzVector(0, 0, 0, 0);
                TLorentzVector* lep3 = new TLorentzVector(0, 0, 0, 0);
                TLorentzVector* lep4 = new TLorentzVector(0, 0, 0, 0);
                TLorentzVector* gph = new TLorentzVector(0, 0, 0, 0);

                // the function returns: lep1, lep3 are the paricles, lep2, lep4 are the anti-particles
                toolbox::fillGenBorn(genPartArr, BOSON_ID, vec, lep1, lep2, lep3, lep4);

                double ptWeight = 1;
                //for (int i = 0; i <= hh_diff->GetNbinsX(); ++i) {
                //    if (vec->Pt() > hh_diff->GetBinLowEdge(i) && vec->Pt() < hh_diff->GetBinLowEdge(i + 1)) {
                //        ptWeight = hh_diff->GetBinContent(i);
                //        break;
                //    }
                //}

                if (charge == 1) {
                    //For W+->e+vu decay, change things up so that lep1 and lep3 are the charged particles
                    TLorentzVector* tmp = lep1;
                    lep1 = lep2;
                    lep2 = tmp;
                    tmp = lep3;
                    lep3 = lep4;
                    lep4 = tmp;
                }

                if (doDressed) {
                    for (Int_t i = 0; i < genPartArr->GetEntriesFast(); i++) {
                        const baconhep::TGenParticle* genloop = (baconhep::TGenParticle*)((*genPartArr)[i]);
                        if (fabs(genloop->pdgId) != 22)
                            continue;
                        gph->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                        if (toolbox::deltaR(gph->Eta(), gph->Phi(), lep3->Eta(), lep3->Phi()) < 0.1) {
                            lep3->operator+=(*gph);
                        }
                    }
                }

                Double_t weight = (gen->weight > 0 ? 1 : -1) * samp->xsecv[ifile] / nWgtSum;
                nEvtsv += weight;
                nEvtsv_pT += weight * ptWeight;

                nEvtsv_QCD[0] += weight * gen->lheweight[1];
                nEvtsv_QCD[1] += weight * gen->lheweight[2];
                nEvtsv_QCD[2] += weight * gen->lheweight[3];
                nEvtsv_QCD[3] += weight * gen->lheweight[4];
                nEvtsv_QCD[4] += weight * gen->lheweight[6];
                nEvtsv_QCD[5] += weight * gen->lheweight[8];

                for (int npdf = 0; npdf < NPDF; npdf++)
                    nEvtsv_PDF[npdf] += weight * gen->lheweight[9 + npdf];

                Bool_t isBarrel = kTRUE;
                if (doborn) {
                    if (lep1->Pt() < PT_CUT)
                        continue;
                    if (fabs(lep1->Eta()) > ETA_CUT)
                        continue;
                    isBarrel = (fabs(lep1->Eta()) < ETA_BARREL) ? kTRUE : kFALSE;
                } else {
                    if (lep3->Pt() < PT_CUT)
                        continue;
                    if (fabs(lep3->Eta()) > ETA_CUT)
                        continue;
                    isBarrel = (fabs(lep3->Eta()) < ETA_BARREL) ? kTRUE : kFALSE;
                }

                nEvtsAfter1Lep += weight;

                double mtgen = 0;
                if (doborn)
                    mtgen = sqrt(2.0 * (lep1->Pt()) * (lep2->Pt()) * (1.0 - cos(toolbox::deltaPhi(lep1->Phi(), lep2->Phi()))));
                else
                    mtgen = sqrt(2.0 * (lep3->Pt()) * (lep4->Pt()) * (1.0 - cos(toolbox::deltaPhi(lep3->Phi(), lep4->Phi()))));

                if (mtgen < 20)
                    continue;

                nEvtsAfterMT += weight;

                nSelv += weight;
                nSelv_pT += weight * ptWeight;
                if (isBarrel)
                    nSelBv += weight;
                else
                    nSelEv += weight;

                nSelv_QCD[0] += weight * gen->lheweight[1];
                nSelv_QCD[1] += weight * gen->lheweight[2];
                nSelv_QCD[2] += weight * gen->lheweight[3];
                nSelv_QCD[3] += weight * gen->lheweight[4];
                nSelv_QCD[4] += weight * gen->lheweight[6];
                nSelv_QCD[5] += weight * gen->lheweight[8];
                for (int npdf = 0; npdf < NPDF; npdf++)
                    nSelv_PDF[npdf] += weight * gen->lheweight[9 + npdf];

                delete vec;
                delete lep1;
                delete lep2;
                delete lep3;
                delete lep4;
                delete gph;
            }
        }

        std::cout << "nselv " << nSelv << "  nevtsv " << nEvtsv << std::endl;

        delete infile;
        infile = 0, eventTree = 0;
        samp = 0;
    }

    // compute acceptances
    accv = nSelv / nEvtsv;
    accErrv = sqrt(accv * (1. - accv) / nEvtsv);
    accBv = nSelBv / nEvtsv;
    accErrBv = sqrt(accBv * (1. - accBv) / nEvtsv);
    accEv = nSelEv / nEvtsv;
    accErrEv = sqrt(accEv * (1. - accEv) / nEvtsv);

    accv_pT = nSelv_pT / nEvtsv_pT;

    // Print full set for efficiency calculations
    char masterOutput[600];
    // just start printing....
    sprintf(masterOutput, "%s/%s.txt", outputDir.Data(), outputName.Data());
    ofstream txtfile;
    txtfile.open(masterOutput);
    txtfile << "acc " << nSelv / nEvtsv << endl;

    for (int j = 0; j < NPDF; ++j)
        txtfile << "pdf" << j << " " << nSelv_PDF[j] / nEvtsv_PDF[j] << endl;
    for (int j = 0; j < NQCD; ++j)
        txtfile << "qcd" << j << " " << nSelv_QCD[j] / nEvtsv_QCD[j] << endl;
    txtfile.close();

    delete gen;

    //--------------------------------------------------------------------------------------------------------------
    // Output
    //==============================================================================================================

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

    cout << "    *** Acceptance ***" << endl;
    cout << "            b: " << setw(12) << nSelBv << " / " << nEvtsv << " = " << accBv << " +/- " << accErrBv << endl;
    cout << "            e: " << setw(12) << nSelEv << " / " << nEvtsv << " = " << accEv << " +/- " << accErrEv << endl;
    cout << "        total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv << endl;
    cout << " with pt: " << setw(12) << accv_pT << endl;
    cout << " pt diff: " << setw(12) << 100 * fabs(accv / accv_pT - 1) << endl;
    cout << endl;

    char txtfname1[300];
    sprintf(txtfname1, "%s/gen.txt", outputDir.Data());
    ofstream txtfile1;
    txtfile1.open(txtfname1);
    txtfile1 << "*" << endl;
    txtfile1 << "* SUMMARY" << endl;
    txtfile1 << "*--------------------------------------------------" << endl;
    if (charge == 0)
        txtfile1 << " W -> mu nu" << endl;
    if (charge == -1)
        txtfile1 << " W- -> mu nu" << endl;
    if (charge == 1)
        txtfile1 << " W+ -> mu nu" << endl;
    txtfile1 << "  pT > " << PT_CUT << endl;
    txtfile1 << "  |eta| < " << ETA_CUT << endl;
    txtfile1 << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
    txtfile1 << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
    txtfile1 << endl;

    txtfile1 << "    *** Acceptance ***" << endl;
    txtfile1 << "            b: " << setw(12) << nSelBv << " / " << nEvtsv << " = " << accBv << " +/- " << accErrBv << endl;
    txtfile1 << "            e: " << setw(12) << nSelEv << " / " << nEvtsv << " = " << accEv << " +/- " << accErrEv << endl;
    txtfile1 << "        total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv << endl;
    txtfile1 << " with pt: " << setw(12) << accv_pT << endl;
    txtfile1 << " pt diff: " << setw(12) << 100 * fabs(accv / accv_pT - 1) << endl;
    txtfile1 << endl;
    txtfile1.close();

    char txtfname2[300];
    sprintf(txtfname2, "%s/acceptance.txt", outputDir.Data());
    ofstream txtfile2;
    txtfile2.open(txtfname2);
    double ndiv = nEvtsv;
    if (charge == 1)
        txtfile2 << "\\PW^{+}\\to\\mu^{+}\\nu" << endl;
    else if (charge == -1)
        txtfile2 << "\\PW^{-}\\to\\mu^{-}\\nu" << endl;
    else
        txtfile2 << "\\PW^{\\pm}\\to\\mu^{\\pm}\\nu" << endl;
    txtfile2 << " Total " << setw(20) << nEvtsv << setw(20) << nEvtsv / ndiv << endl;
    txtfile2 << " After_lep1_cut " << setw(20) << nEvtsAfter1Lep << setw(20) << nEvtsAfter1Lep / ndiv << endl;
    txtfile2 << " After_MT_cut " << setw(20) << nEvtsAfterMT << setw(20) << nEvtsAfterMT / ndiv << endl;
    txtfile2 << endl;
    txtfile2.close();

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;

    gBenchmark->Show("computeAccGenWm");
}
