//================================================================================================
//
// Compute Z->ll acceptance at generator level
//
//  * outputs results summary text file
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TBenchmark.h>   // class to track macro running statistics
#include <TClonesArray.h> // ROOT array class
#include <TFile.h>        // file handle class
#include <TH1D.h>         // histogram class
#include <TLorentzVector.h>
#include <TROOT.h>   // access to gROOT, entry point to ROOT system
#include <TSystem.h> // interface to OS
#include <TTree.h>   // class to access ntuples
#include <fstream>   // functions for file I/O
#include <iomanip>   // functions to format standard I/O
#include <iostream>  // standard I/O
#include <sstream>   // class for parsing strings
#include <string>    // C++ string class
#include <vector>    // STL vector class

// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "MitEwk13TeV/Utils/CSample.hh"   // helper class to handle samples
#include "MitEwk13TeV/Utils/ConfParse.hh" // input conf file parser
#include "MitEwk13TeV/Utils/MyTools.hh"
#endif

//=== MAIN MACRO =================================================================================================

void computeAccGenZll(const TString conf,       // input file
                      const TString outputDir,  // output directory
                      const TString outputName, // output filename
                      const bool doDressed = 0, // do dressed
                      const bool doborn = 0,    // do born
                      const bool doMuon = 1,    // do muon
                      const float frac = 0.1    // fraction of events for calculation
)
{
    gBenchmark->Start("computeAccGenZll");

    //--------------------------------------------------------------------------------------------------------------
    // Settings
    //==============================================================================================================
    // bool doPtWeights = true;
    const Double_t MASS_LOW = 60;
    const Double_t MASS_HIGH = 120;
    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;
    const Double_t ETA_BARREL = 1.2;
    const Double_t ETA_ENDCAP = 1.2;

    const Int_t BOSON_ID = 23;
    Int_t LEPTON_ID = 13;
    if (!doMuon)
        LEPTON_ID = 11;

    const Int_t NPDF = 100;
    const Int_t NQCD = 6;

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    vector<TString> snamev;    // sample name (for output files)
    vector<CSample *> samplev; // data/MC samples

    //
    // parse .conf file
    //
    confParse(conf, snamev, samplev);

    // Create output directory
    gSystem->mkdir(outputDir, kTRUE);

    // Data structures to store info from TTrees
    baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
    TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");

    TFile *infile = 0;
    TTree *eventTree = 0;

    // Variables to store acceptances and uncertainties
    Double_t nEvtsv = 0;
    Double_t nEvts_noWeight = 0;
    Double_t nEvtsBeforeMass = 0, nEntries = 0, nEvtsAfterMass = 0, nEvtsAfter1Lep = 0, nEvtsAfter2Lep = 0;
    Double_t nSelv = 0, accv = 0, accErrv = 0;
    Double_t nSelBBv = 0, accBBv = 0, accErrBBv = 0;
    Double_t nSelBEv = 0, accBEv = 0, accErrBEv = 0;
    Double_t nSelEEv = 0, accEEv = 0, accErrEEv = 0;

    Double_t nEvtsv_pT = 0, nSelv_pT = 0, accv_pT = 0;

    vector<Double_t> nEvtsv_QCD(NQCD, 0), nSelv_QCD(NQCD, 0);
    vector<Double_t> nEvtsv_PDF(NPDF, 0), nSelv_PDF(NPDF, 0);

    TString sqrts = "13TeV";
    if (conf.Contains("5"))
        sqrts = "5TeV";
    // TFile* rf = new TFile("/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/SignalExtraction/Z_pT/zPt_Normal" + sqrts + ".root");
    // TH1D* hh_diff = (TH1D*)rf->Get("hZptRatio");

    //
    // loop through files
    //
    CSample *samp = samplev[0];
    const UInt_t nfiles = samp->fnamev.size();

    for (UInt_t ifile = 0; ifile < nfiles; ifile++)
    {

        // Read input file and get the TTrees
        cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... ";
        cout.flush();

        infile = TFile::Open(samp->fnamev[ifile]);
        assert(infile);

        eventTree = (TTree *)infile->Get("Events");
        assert(eventTree);
        eventTree->SetBranchAddress("GenEvtInfo", &gen);
        TBranch *genBr = eventTree->GetBranch("GenEvtInfo");
        eventTree->SetBranchAddress("GenParticle", &genPartArr);
        TBranch *partBr = eventTree->GetBranch("GenParticle");

        //
        // loop over events
        //
        double nWgtSum = 0., nAbsSum = 0; // total number of events after reweighting

        // loop over the events first, to get the positive and negative frations of events,
        // used for scaling later
        std::cout << "Process events quickly first time for counting" << std::endl;
        for (UInt_t ientry = 0; ientry < (uint)(frac * eventTree->GetEntries()); ientry++)
        {
            if (ientry % 100000 == 0)
                cout << "Processing event " << ientry << ". " << (double)ientry / (double)eventTree->GetEntries() * 100 << " percent done with this file." << endl;
            genBr->GetEntry(ientry);
            nAbsSum += 1.;
            nWgtSum += (gen->weight > 0) ? 1 : -1;
        }
        std::cout << "Finished first loop. Total events " << nAbsSum << " after negative weight subtraction " << nWgtSum / nAbsSum << std::endl;

        for (UInt_t ientry = 0; ientry < (uint)(frac * eventTree->GetEntries()); ientry++)
        {
            if (ientry % 100000 == 0)
                cout << "Processing event " << ientry << ". " << (double)ientry / (double)eventTree->GetEntries() * 100 << " percent done with this file." << endl;

            genBr->GetEntry(ientry);
            genPartArr->Clear();
            partBr->GetEntry(ientry);

            Int_t lepq1 = -99, lepq2 = -99;

            if (fabs(toolbox::flavor(genPartArr, BOSON_ID)) != LEPTON_ID)
                continue;
            nEntries += 1;

            TLorentzVector *vec = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector *lep1 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector *lep2 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector *lep3 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector *lep4 = new TLorentzVector(0, 0, 0, 0);
            TLorentzVector *gph = new TLorentzVector(0, 0, 0, 0);

            toolbox::fillGenBorn(genPartArr, BOSON_ID, vec, lep1, lep2, lep3, lep4);

            double ptWeight = 1;
            // for (int i = 0; i <= hh_diff->GetNbinsX(); ++i) {
            //     if (vec->Pt() > hh_diff->GetBinLowEdge(i) && vec->Pt() < hh_diff->GetBinLowEdge(i + 1)) {
            //         ptWeight = hh_diff->GetBinContent(i);
            //         break;
            //     }
            // }

            if (doDressed)
            {
                for (Int_t i = 0; i < genPartArr->GetEntriesFast(); i++)
                {
                    const baconhep::TGenParticle *genloop = (baconhep::TGenParticle *)((*genPartArr)[i]);
                    if (fabs(genloop->pdgId) != 22)
                        continue;
                    gph->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                    if (toolbox::deltaR(gph->Eta(), gph->Phi(), lep3->Eta(), lep3->Phi()) < 0.1)
                        lep3->operator+=(*gph);
                    if (toolbox::deltaR(gph->Eta(), gph->Phi(), lep4->Eta(), lep4->Phi()) < 0.1)
                        lep4->operator+=(*gph);
                }
            }

            Double_t weight = (gen->weight > 0 ? 1 : -1) * samp->xsecv[ifile] / nWgtSum;

            nEvtsBeforeMass += weight;

            // double genmass = vec->M();
            TLorentzVector dilep;
            if (doborn)
                dilep = (*lep1) + (*lep2);
            else
                dilep = (*lep3) + (*lep4);

            if (dilep.M() < MASS_LOW || dilep.M() > MASS_HIGH)
                continue;

            nEvtsAfterMass += weight;

            nEvts_noWeight += 1;
            nEvtsv += weight;
            nEvtsv_pT += weight * ptWeight;

            // -------------------------------------------------
            // I'm not sure which indexing is correct for Z's
            // -------------------------------------------------
            nEvtsv_QCD[0] += weight * gen->lheweight[0];
            nEvtsv_QCD[1] += weight * gen->lheweight[1];
            nEvtsv_QCD[2] += weight * gen->lheweight[2];
            nEvtsv_QCD[3] += weight * gen->lheweight[3];
            nEvtsv_QCD[4] += weight * gen->lheweight[5];
            nEvtsv_QCD[5] += weight * gen->lheweight[7];

            for (int npdf = 0; npdf < NPDF; npdf++)
                nEvtsv_PDF[npdf] += weight * gen->lheweight[8 + npdf];

            Bool_t isB1, isB2;
            if (doborn)
            {
                if (lep1->Pt() < PT_CUT)
                    continue;
                if (fabs(lep1->Eta()) > ETA_CUT)
                    continue;
                isB1 = (fabs(lep1->Eta()) < ETA_BARREL) ? kTRUE : kFALSE;

                nEvtsAfter1Lep += weight;

                if (lep2->Pt() < PT_CUT)
                    continue;
                if (fabs(lep2->Eta()) > ETA_CUT)
                    continue;
                isB2 = (fabs(lep2->Eta()) < ETA_BARREL) ? kTRUE : kFALSE;

                nEvtsAfter2Lep += weight;
            }
            else
            {
                if (lep3->Pt() < PT_CUT)
                    continue;
                if (fabs(lep3->Eta()) > ETA_CUT)
                    continue;
                isB1 = (fabs(lep3->Eta()) < ETA_BARREL) ? kTRUE : kFALSE;

                nEvtsAfter1Lep += weight;

                if (lep4->Pt() < PT_CUT)
                    continue;
                if (fabs(lep4->Eta()) > ETA_CUT)
                    continue;
                isB2 = (fabs(lep4->Eta()) < ETA_BARREL) ? kTRUE : kFALSE;

                nEvtsAfter2Lep += weight;
            }

            nSelv += weight;
            nSelv_pT += weight * ptWeight;
            if (isB1 && isB2)
                nSelBBv += weight;
            else if (!isB1 && !isB2)
                nSelEEv += weight;
            else
                nSelBEv += weight;

            nSelv_QCD[0] += weight * gen->lheweight[0];
            nSelv_QCD[1] += weight * gen->lheweight[1];
            nSelv_QCD[2] += weight * gen->lheweight[2];
            nSelv_QCD[3] += weight * gen->lheweight[3];
            nSelv_QCD[4] += weight * gen->lheweight[5];
            nSelv_QCD[5] += weight * gen->lheweight[7];
            for (int npdf = 0; npdf < NPDF; npdf++)
                nSelv_PDF[npdf] += weight * gen->lheweight[8 + npdf];

            delete vec;
            delete lep1;
            delete lep2;
            delete lep3;
            delete lep4;
            delete gph;
        }

        // compute acceptances
        // unc = sqrt((acc * (1. - acc)) / N);
        accv = (nSelv / nEvtsv);
        accErrv = (sqrt(accv * (1. - accv) / nEvts_noWeight));
        accBBv = (nSelBBv / nEvtsv);
        accErrBBv = (sqrt(accBBv * (1. - accBBv) / nEvts_noWeight));
        accBEv = (nSelBEv / nEvtsv);
        accErrBEv = (sqrt(accBEv * (1. - accBEv) / nEvts_noWeight));
        accEEv = (nSelEEv / nEvtsv);
        accErrEEv = (sqrt(accEEv * (1. - accEEv) / nEvts_noWeight));

        accv_pT = (nSelv_pT / nEvtsv_pT);

        std::cout << "nselv " << nSelv << "  nevtsv " << nEvtsv << std::endl;

        delete infile;
        infile = 0, eventTree = 0;
    }

    delete gen;

    char masterOutput[600];
    sprintf(masterOutput, "%s/%s.txt", outputDir.Data(), outputName.Data());
    ofstream txtfile;
    txtfile.open(masterOutput);
    txtfile << "acc " << nSelv / nEvtsv << endl;
    for (int j = 0; j < NPDF; ++j)
        txtfile << "pdf" << j << " " << nSelv_PDF[j] / nEvtsv_PDF[j] << endl;
    for (int j = 0; j < NQCD; ++j)
        txtfile << "qcd" << j << " " << nSelv_QCD[j] / nEvtsv_QCD[j] << endl;
    txtfile.close();

    //--------------------------------------------------------------------------------------------------------------
    // Output
    //==============================================================================================================
    cout << "*" << endl;
    cout << "* SUMMARY" << endl;
    cout << "*--------------------------------------------------" << endl;
    if (doMuon)
        cout << " Z -> mu mu" << endl;
    else
        cout << " Z -> e e" << endl;
    cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
    cout << "  pT > " << PT_CUT << endl;
    cout << "  |eta| < " << ETA_CUT << endl;
    cout << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
    cout << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
    cout << endl;

    cout << " total entries " << nEntries << endl;
    cout << " total events after GEN mass cut " << nEvtsAfterMass << endl;
    cout << " total events before the mass cut " << nEvtsBeforeMass << endl;
    cout << " total events after the gen mass cut " << nEvtsv << endl;
    cout << " total events after the gen lepton pt and eta cut " << nSelv << endl;
    cout << "    *** Acceptance ***" << endl;
    cout << "     b-b: " << setw(12) << nSelBBv << " / " << nEvtsv << " = " << accBBv << " +/- " << accErrBBv << endl;
    cout << "     b-e: " << setw(12) << nSelBEv << " / " << nEvtsv << " = " << accBEv << " +/- " << accErrBEv << endl;
    cout << "     e-e: " << setw(12) << nSelEEv << " / " << nEvtsv << " = " << accEEv << " +/- " << accErrEEv << endl;
    cout << "   total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv << endl;
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
    if (doMuon)
        txtfile1 << " Z -> mu mu" << endl;
    else
        txtfile1 << " Z -> e e" << endl;
    txtfile1 << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
    txtfile1 << "  pT > " << PT_CUT << endl;
    txtfile1 << "  |eta| < " << ETA_CUT << endl;
    txtfile1 << "  Barrel definition: |eta| < " << ETA_BARREL << endl;
    txtfile1 << "  Endcap definition: |eta| > " << ETA_ENDCAP << endl;
    txtfile1 << endl;

    txtfile1 << "    *** Acceptance ***" << endl;
    txtfile1 << "     b-b: " << setw(12) << nSelBBv << " / " << nEvtsv << " = " << accBBv << " +/- " << accErrBBv << endl;
    txtfile1 << "     b-e: " << setw(12) << nSelBEv << " / " << nEvtsv << " = " << accBEv << " +/- " << accErrBEv << endl;
    txtfile1 << "     b-e: " << setw(12) << nSelEEv << " / " << nEvtsv << " = " << accEEv << " +/- " << accErrEEv << endl;
    txtfile1 << "   total: " << setw(12) << nSelv << " / " << nEvtsv << " = " << accv << " +/- " << accErrv << endl;
    txtfile1 << " with pt: " << setw(12) << accv_pT << endl;
    txtfile1 << " pt diff: " << setw(12) << 100 * fabs(accv / accv_pT - 1) << endl;
    txtfile1 << endl;
    txtfile1.close();

    char txtfname2[300];
    sprintf(txtfname2, "%s/acceptance.txt", outputDir.Data());
    ofstream txtfile2;
    txtfile2.open(txtfname2);
    double ndiv = nEvtsAfterMass;
    if (doMuon)
        txtfile2 << " \\PZ\\to\\mu^{+}\\mu^{-}" << endl;
    else
        txtfile2 << " \\PZ\\to e^{+}e^{-}" << endl;
    txtfile2 << " Total " << setw(20) << nEvtsAfterMass << setw(20) << nEvtsAfterMass / ndiv << endl;
    txtfile2 << " After_lep1_cut " << setw(20) << nEvtsAfter1Lep << setw(20) << nEvtsAfter1Lep / ndiv << endl;
    txtfile2 << " After_lep2_cut " << setw(20) << nEvtsAfter2Lep << setw(20) << nEvtsAfter2Lep / ndiv << endl;
    txtfile2 << endl;
    txtfile2.close();

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;

    gBenchmark->Show("computeAccGenZll");
}
