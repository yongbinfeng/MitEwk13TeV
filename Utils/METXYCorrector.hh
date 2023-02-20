#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TTree.h"
#include <sstream>
#include <string>
#include <vector>

#include <math.h>
#include <stdio.h>
//
// ** apply MET XY corrections **
//
// code is modified based on https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18andUL16.h
//

using namespace std;

class METXYCorrector
{

public:
    METXYCorrector(std::string iName, std::string iFName);
    void loadXYCorrection(string iNameFile);
    void CorrectMETXY(float &pfmet, float &pfmetphi, bool isData);

private:
    float mean_met_x_data;
    float mean_met_y_data;
    float mean_met_x_mc;
    float mean_met_y_mc;
};

//-----------------------------------------------------------------------------------------------------------------------------------------
// constructors, but some of the parts are pointless lol

METXYCorrector::METXYCorrector(std::string iName, std::string iFName)
{
    loadXYCorrection(iFName);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
void METXYCorrector::loadXYCorrection(std::string iFName)
{
    TFile *lFile = new TFile(iFName.c_str());
    mean_met_x_data = ((TH1D *)lFile->Get("hMET_x_Data"))->GetMean();
    mean_met_y_data = ((TH1D *)lFile->Get("hMET_y_Data"))->GetMean();
    mean_met_x_mc = ((TH1D *)lFile->Get("hMET_x_MC"))->GetMean();
    mean_met_y_mc = ((TH1D *)lFile->Get("hMET_y_MC"))->GetMean();
    lFile->Delete();
    std::cout << "load XY correction from " << iFName << std::endl;
    std::cout << "mean_met_x_data " << mean_met_x_data << " mean_met_y_data " << mean_met_y_data << std::endl;
    std::cout << "mean_met_x_mc " << mean_met_x_mc << " mean_met_y_mc " << mean_met_y_mc << std::endl;
    std::cout << " ========= " << std::endl;
}

void METXYCorrector::CorrectMETXY(float &pfmet, float &pfmetphi, bool isData)
{
    float pfmetx = pfmet * TMath::Cos(pfmetphi);
    float pfmety = pfmet * TMath::Sin(pfmetphi);

    // std::cout << "pfmetx " << pfmetx << " pfmety " << pfmety << " nPV " << nPV;
    if (isData)
    {
        pfmetx -= mean_met_x_data;
        pfmety -= mean_met_y_data;
    }
    else
    {
        pfmetx -= mean_met_x_mc;
        pfmety -= mean_met_y_mc;
    }
    // std::cout << " after corr pfmetx " << pfmetx << " pfmety " << pfmety << std::endl;

    pfmet = TMath::Sqrt(pfmetx * pfmetx + pfmety * pfmety);
    if (pfmetx == 0 && pfmety > 0)
        pfmetphi = TMath::Pi();
    else if (pfmetx == 0 && pfmety < 0)
        pfmetphi = -TMath::Pi();
    else if (pfmetx > 0)
        pfmetphi = TMath::ATan(pfmety / pfmetx);
    else if (pfmetx < 0 && pfmety > 0)
        pfmetphi = TMath::ATan(pfmety / pfmetx) + TMath::Pi();
    else if (pfmetx < 0 && pfmety < 0)
        pfmetphi = TMath::ATan(pfmety / pfmetx) - TMath::Pi();
    else
        pfmetphi = 0;
}