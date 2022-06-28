#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../../Utils/CPlot.hh" // helper class for plots
#include "../../Utils/MitStyleRemix.hh" // style settings for drawing
#include <TBenchmark.h> // class to track macro running statistics
#include <TCanvas.h> // class for drawing
#include <TEfficiency.h> // class to handle efficiency calculations
#include <TFile.h> // file handle class
#include <TH1D.h> // 1D histograms
#include <TH2D.h> // 2D histograms
#include <TLorentzVector.h> // class for 4-vector calculations
#include <TROOT.h> // access to gROOT, entry point to ROOT system
#include <TStyle.h> // class to handle ROOT plotting style
#include <TSystem.h> // interface to OS
#include <TTree.h> // class to access ntuples
#include <fstream> // functions for file I/O
#include <iomanip> // functions to format standard I/O
#include <iostream> // standard I/O
#include <sstream> // class for parsing strings
#include <string> // C++ string class
#include <vector> // STL vector class
#endif

// RooFit headers
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooMCStudy.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#define BIN_SIZE_PASS 2
#define BIN_SIZE_FAIL 2

void toyGenAndPull_TagPt(const TString tmpDir, // template for generating pseudodata (should be the alternate model)
        const TString fitDir, // templates for running one set of fits (should be the default model)
        const TString binName,
        const TString outputDir,
        const TString outputName,
        const Int_t nPsExp,
        const TString tmpDir2)
{
    std::cout << "fit Dir " << fitDir << std::endl;
    std::cout << "tmp Dir " << tmpDir << std::endl;
    std::cout << "tmp2 Dir " << tmpDir2 << std::endl;

    gSystem->Load("/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Utils/RooCMSShape_cc.so");

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    gSystem->mkdir(outputDir, kTRUE);

    // RETRIEVE TRUTH PDFS
    TFile* ftmp = new TFile(tmpDir + binName + ".root");
    TFile* ftmp2 = new TFile(tmpDir2 + binName + ".root");
    TFile* ffit = new TFile(fitDir + binName + ".root");
    RooWorkspace* wtmp = (RooWorkspace*)ftmp->Get("w");
    RooWorkspace* wtmp2 = (RooWorkspace*)ftmp2->Get("w");
    RooWorkspace* wfit = (RooWorkspace*)ffit->Get("w");

    RooCategory sample("sample", "");
    sample.defineType("Pass", 1);
    sample.defineType("Fail", 2);

    RooRealVar* m = wfit->var("m");

    Double_t effBdata = wtmp->var("eff")->getVal();
    Double_t errBdata = (fabs(wtmp->var("eff")->getErrorLo()) + fabs(wtmp->var("eff")->getErrorHi())) * 0.5;
    Double_t effAdata = wfit->var("eff")->getVal();
    Double_t errAdata = (fabs(wfit->var("eff")->getErrorLo()) + fabs(wfit->var("eff")->getErrorHi())) * 0.5;
   
    // alternate model to generate pseudo data
    RooAbsPdf* modelFail = (RooAbsPdf*)wtmp->pdf("modelFail")->Clone("modelFailGeneration");
    RooAbsPdf* modelPass = (RooAbsPdf*)wtmp->pdf("modelPass")->Clone("modelPassGeneration");
    RooSimultaneous templatePdf("templatePdf", "templatePdf", sample);
    templatePdf.addPdf(*modelPass, "Pass");
    templatePdf.addPdf(*modelFail, "Fail");

    RooAbsPdf* modelFail2 = (RooAbsPdf*)wtmp2->pdf("modelFail")->Clone("modelFailGeneration2");
    RooAbsPdf* modelPass2 = (RooAbsPdf*)wtmp2->pdf("modelPass")->Clone("modelPassGeneration2");
    RooSimultaneous templatePdf2("templatePdf2", "templatePdf2", sample);
    templatePdf2.addPdf(*modelPass2, "Pass");
    templatePdf2.addPdf(*modelFail2, "Fail");

    char outHistName[20];
    sprintf(outHistName, "diffHist_%s", binName.Data());
    char outEffAHistName[20];
    sprintf(outEffAHistName, "effAHist_%s", binName.Data());
    char outEffBHistName[20];
    sprintf(outEffBHistName, "effBHist_%s", binName.Data());
    TH1D* hdiff = new TH1D(outHistName, outHistName, 100, -0.1, 0.1);
    TH1D* heffA = new TH1D(outEffAHistName, outEffAHistName, 100, -5.0, 5.0);
    TH1D* heffB = new TH1D(outEffBHistName, outEffBHistName, 100, -5.0, 5.0);

    TTree* intree = (TTree*)ftmp->Get("Bin");
    UInt_t nEvents;
    intree->SetBranchAddress("nEvents", &nEvents);
    intree->GetEntry(0);

    TTree* intree2 = (TTree*)ftmp2->Get("Bin");
    UInt_t nEvents2;
    intree2->SetBranchAddress("nEvents2", &nEvents2);
    intree2->GetEntry(0);

    // fit model from wfit
    RooAbsPdf* fitmodelFail = (RooAbsPdf*)wfit->pdf("modelFail")->Clone("modelFailFitA");
    RooAbsPdf* fitmodelPass = (RooAbsPdf*)wfit->pdf("modelPass")->Clone("modelPassFitA");
    RooSimultaneous fitPdf("fitPdf", "fitPdf", sample);
    fitPdf.addPdf(*fitmodelPass, "Pass");
    fitPdf.addPdf(*fitmodelFail, "Fail");

    // alternative model from wtemp
    RooWorkspace* wfit2 = (RooWorkspace*)wtmp->Clone("wfit2");
    RooAbsPdf* fitmodel2Fail = (RooAbsPdf*)wfit2->pdf("modelFail")->Clone("modelFailFitB");
    RooAbsPdf* fitmodel2Pass = (RooAbsPdf*)wfit2->pdf("modelPass")->Clone("modelPassFitB");
    // alternate model to generate pseudo data
    RooSimultaneous fit2Pdf("fit2Pdf", "fit2Pdf", sample);
    fit2Pdf.addPdf(*fitmodel2Pass, "Pass");
    fit2Pdf.addPdf(*fitmodel2Fail, "Fail");
    
    RooDataSet *dataCombined = 0;
    RooDataSet *dataCombined2 = 0;
    RooDataSet *dataFit = 0;

    // RooFitResult *fitResult=0;
    std::cout << "templatePdf print " << std::endl;
    templatePdf.Print(); // alternate model to generate data
    std::cout << "templatePdf2 print " << std::endl;
    templatePdf2.Print();

    std::cout << "fitPdf print " << std::endl;
    fitPdf.Print();  // default model to fit the data
    std::cout << "fit2Pdf print " << std::endl;
    fit2Pdf.Print();

    for(int i = 0; i < nPsExp; ++i){
        dataCombined = templatePdf.generate(RooArgList(*m,sample),nEvents,Extended());
        dataCombined2 = templatePdf2.generate(RooArgList(*m, sample), nEvents2, Extended());

        // fit with the alternate model
        RooFitResult *fitResult2 = fit2Pdf.fitTo(*dataCombined,
                RooFit::PrintEvalErrors(-1),
                RooFit::Extended(),
                RooFit::Strategy(2),
                // RooFit::Minos(RooArgSet(*eff)),
                RooFit::Save());

        RooRealVar* eff2 = (RooRealVar*)fitResult2->floatParsFinal().find("eff");
        RooRealVar* eff20 = (RooRealVar*)fitResult2->floatParsInit().find("eff");

        if((fabs(eff2->getErrorLo())<5e-5) || (eff2->getErrorHi()<5e-5) || fitResult2->status()!=0){
            fitResult2 = fit2Pdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());
        }

        dataFit = (RooDataSet*)dataCombined->Clone("dataForFit");
        std::cout << "data combined " << std::endl;
        dataFit->Print("v");
        std::cout << std::endl << std::endl << " data after combining " << std::endl;
        dataFit->append(*dataCombined2);
        dataFit->Print("v");

        RooMsgService::instance().setSilentMode(kTRUE);

        // fit with the default model
        RooFitResult *fitResult = fitPdf.fitTo(*dataFit,
                RooFit::PrintEvalErrors(-1),
                RooFit::Extended(),
                RooFit::Strategy(2),
                // RooFit::Minos(RooArgSet(*eff)),
                RooFit::Save());

        RooRealVar* eff = (RooRealVar*)fitResult->floatParsFinal().find("eff");
        RooRealVar* eff0 = (RooRealVar*)fitResult->floatParsInit().find("eff");

        if((fabs(eff->getErrorLo())<5e-5) || (eff->getErrorHi()<5e-5) || fitResult->status()!=0){
            fitResult = fitPdf.fitTo(*dataFit, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());
        }

        std::cout << "RooFit Status " << fitResult->status() << " alternative " << fitResult2->status() << std::endl;
        std::cout <<"EFFICIENCY A|B " << eff->getVal()<<" A|data "<< effAdata << " B|B " <<  eff2->getVal() << " B|data " <<  effBdata <<  std::endl;
        hdiff->Fill(eff->getVal() - eff2->getVal());
        heffA->Fill((eff->getVal() - effAdata) / errAdata);
        heffB->Fill((eff2->getVal() - effBdata) / errBdata);
        std::cout << "pull eff A " << (eff->getVal() - effAdata) / errAdata << " pull eff B " << (eff2->getVal() - effBdata) / errBdata << std::endl;
        std::cout << "count " << templatePdf.getVal() << std::endl;

    }
    // comment for now

    // std::cout << "blah" << std::endl;

    TF1 *f = new TF1("f","gaus",-0.1, 0.1);
    f->SetLineColor(2);
    f->SetParameters(nPsExp,0.0,0.01);
    hdiff->Fit("f");
    float meanhist = hdiff->GetMean();
    // float meangaus = abs(f->GetParameter(1));
    float meangaus = f->GetParameter(1);
    // float mean = meangaus < meanhist ? meangaus : meanhist;
    // // float mean = meantemp > 0? meantemp : -1.*meantemp;
    float mean = meangaus;
    float sigma = f->GetParameter(2);

    std::cout << "mean hist " << meanhist << " mean gaus " << mean << std::endl;

    char outputfile[100];
    sprintf(outputfile,"%s/%s.txt",outputDir.Data(), outputName.Data());
    ofstream meanfile;
    meanfile.open(outputfile);
    // gaus mean
    meanfile << mean <<endl;
    // hist mean
    meanfile << meanhist << endl;
    meanfile.close();

    sprintf(outputfile, "%s/effDiff_%s.png",outputDir.Data(), outputName.Data());
    TCanvas *c = new TCanvas("c","c");
    //gStyle->SetOptFit();
    hdiff->GetXaxis()->SetTitle("#varepsilon_{A|B} - #varepsilon_{B|B}");
    hdiff->Draw();
    c->SaveAs(outputfile);
    c->Clear();

    sprintf(outputfile, "%s/effA_%s.png",outputDir.Data(), outputName.Data());
    TCanvas *ceffA = new TCanvas("ceffA", "ceffA");
    heffA->SetLineColor(3);
    heffA->GetXaxis()->SetTitle("#varepsilon Pull: (#varepsilon_{X|B} - #varepsilon_{X|data}) / statUnc");
    heffA->SetMaximum(1.3*max(heffA->GetMaximum(), heffB->GetMaximum()));
    heffA->Draw();
    heffB->SetLineColor(4);
    heffB->Draw("same");
    ceffA->SaveAs(outputfile);
    ceffA->Clear();
}
