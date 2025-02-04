#ifndef App_Eff_SF
#define App_Eff_SF

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
#include <fstream>
#include <fstream>  // functions for file I/O
#include <iomanip>  // functions to format standard I/O
#include <iostream> // standard I/O
#include <sstream>
#include <sstream> // class for parsing strings
#include <stdexcept>
#include <string> // C++ string class
#include <vector> // STL vector class
// helper class to handle efficiency tables
#include "MitEwk13TeV/Utils/CEffUser2D.hh"

struct basicEff
{
    // data pos
    CEffUser2D dataPos;
    // data neg
    CEffUser2D dataNeg;
    // mc pos
    CEffUser2D mcPos;
    // mc neg
    CEffUser2D mcNeg;
};

struct effUnc
{
    // fsr pos
    CEffUser2D fsrPos;
    // fsr neg
    CEffUser2D fsrNeg;
    // mc pos
    CEffUser2D mcPos;
    // mc neg
    CEffUser2D mcNeg;
    // bkg pos
    CEffUser2D bkgPos;
    // bkg neg
    CEffUser2D bkgNeg;
    // tag pt pos
    CEffUser2D tagPos;
    // tag pt neg
    CEffUser2D tagNeg;
};

// Make everything public for now while I test this
class AppEffSF
{
public:
    AppEffSF() {}
    ~AppEffSF() {}
    AppEffSF(TString dirName)
    {
        this->dirName = dirName;
        isMuon = false;
    }

    void SetBaseDir(TString dirName)
    {
        this->dirName = dirName;
    }

    // some functions to load the individual correction locations
    void loadHLT(TString dir, TString q1, TString q2)
    {
        loadEfficiencyFile(hlt, dir, q1, q2);
        return;
    }

    void loadSta(TString dir, TString q1, TString q2)
    {
        loadEfficiencyFile(sta, dir, q1, q2);
        isMuon = true;
        return;
    }

    void loadSel(TString dir, TString q1, TString q2)
    {
        loadEfficiencyFile(sel, dir, q1, q2);
        return;
    }

    void loadUncSel(TString dir)
    {
        cout << "Loading Selection Uncertainties " << endl;
        loadUncertaintyFile(u_sel, dir);
        return;
    }

    void loadUncSta(TString dir)
    {
        loadUncertaintyFile(u_sta, dir);
        return;
    }

    // do i also need for hlt?
    vector<double> getUncSta(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        return computeAllUncertainties(sta, u_sta, l1, q1, l2, q2);
    }

    vector<double> getUncSel(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        return computeAllUncertainties(sel, u_sel, l1, q1, l2, q2);
    }

    // return a vector with each of the uncertainties propagated to final acc values
    vector<double> computeAllUncertainties(basicEff &effs, effUnc &uncs, TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        vector<double> res;
        res.push_back(computeWithUnc(effs, uncs.fsrPos, uncs.fsrNeg, l1, q1, l2, q2));

        res.push_back(computeWithUnc(effs, uncs.mcPos, uncs.mcNeg, l1, q1, l2, q2));
        res.push_back(computeWithUnc(effs, uncs.bkgPos, uncs.bkgNeg, l1, q1, l2, q2));
        res.push_back(computeWithUnc(effs, uncs.tagPos, uncs.tagNeg, l1, q1, l2, q2));

        return res;
    }

    double computeWithUnc(basicEff &effs, CEffUser2D &varPos, CEffUser2D &varNeg, TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        double effdata = 1, effmc = 1;
        if (q1 > 0)
        {
            effmc *= effs.mcPos.getEff(l1->Eta(), l1->Pt());
            effdata *= effs.dataPos.getEff(l1->Eta(), l1->Pt()) * varPos.getEff(l1->Eta(), l1->Pt());
        }
        else
        {
            effmc *= effs.mcNeg.getEff(l1->Eta(), l1->Pt());
            effdata *= effs.dataNeg.getEff(l1->Eta(), l1->Pt()) * varNeg.getEff(l1->Eta(), l1->Pt());
        }
        if (!l2)
        {
            return effdata / effmc;
        }
        if (q2 > 0)
        {
            effmc *= effs.mcPos.getEff(l2->Eta(), l2->Pt());
            effdata *= effs.dataPos.getEff(l2->Eta(), l2->Pt()) * varPos.getEff(l2->Eta(), l2->Pt());
        }
        else
        {
            effmc *= effs.mcNeg.getEff(l2->Eta(), l2->Pt());
            effdata *= effs.dataNeg.getEff(l2->Eta(), l2->Pt()) * varNeg.getEff(l2->Eta(), l2->Pt());
        }
        return effdata / effmc;
    }

    // compute the corrections from mc to data: eff_data / eff_mc
    double fullCorrections(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        // if there is a second lepton we need to calculate twice
        double corr = 1;
        corr *= computeHLTSF(l1, q1, l2, q2);
        corr *= computeSelSF(l1, q1, l2, q2);
        if (isMuon)
        {
            corr *= computeStaSF(l1, q1, l2, q2);
        }

        return corr;
    }

    double efficiencyOnly(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0, bool isData = 1)
    {
        // compute the efficiencies on data or on mc
        double eff = 1;
        eff *= effHLT(l1, q1, l2, q2, isData);
        eff *= effSel(l1, q1, l2, q2, isData);
        if (isMuon)
        {
            eff *= effSta(l1, q1, l2, q2, isData);
        }
        return eff;
    }

    double statHLTdilep(TH2D *pos, TH2D *neg, double wgt, TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        double stat1 = statUnc(hlt, l1, q1, pos, neg, wgt);
        double stat2 = statUnc(hlt, l2, q2, pos, neg, wgt);
        double eff1 = 1 - computeHLTSF(l1, q1);
        double eff2 = 1 - computeHLTSF(l2, q2);
        cout << stat1 << " " << stat2 << " " << eff1 << " " << eff2 << endl;
        return sqrt(stat1 * stat1 * eff1 * eff1 + stat2 * stat2 * eff2 * eff2);
    }

    double statUncHLT(TLorentzVector *l1, int q1, TH2D *pos, TH2D *neg, double wgt, bool scale = false)
    {
        return statUnc(hlt, l1, q1, pos, neg, wgt, scale);
    }
    double statUncSel(TLorentzVector *l1, int q1, TH2D *pos, TH2D *neg, double wgt, bool scale = false)
    {
        return statUnc(sel, l1, q1, pos, neg, wgt, scale);
    }
    double statUncSta(TLorentzVector *l1, int q1, TH2D *pos, TH2D *neg, double wgt, bool scale = false)
    {
        return statUnc(sta, l1, q1, pos, neg, wgt, scale);
    }

    // private:

    // structs i'm storying the corrections in
    basicEff hlt;
    basicEff sta;
    basicEff sel;

    effUnc u_hlt;
    effUnc u_sta;
    effUnc u_sel;

    TString dirName; // name of primary directory
    bool isMuon;

    void loadEfficiencyFile(basicEff &effs, TString dir, TString q1, TString q2)
    {

        TString dataPosName = dirName + "Data/" + dir + "/" + q1 + "/eff.root";
        TString dataNegName = dirName + "Data/" + dir + "/" + q2 + "/eff.root";
        TString mcPosName = dirName + "MC/" + dir + "/" + q1 + "/eff.root";
        TString mcNegName = dirName + "MC/" + dir + "/" + q2 + "/eff.root";
        cout << "attempting to load " << dataPosName << endl;

        TFile *dp = new TFile(dataPosName);
        effs.dataPos.loadEff((TH2D *)dp->Get("hEffEtaPt"), (TH2D *)dp->Get("hErrlEtaPt"), (TH2D *)dp->Get("hErrhEtaPt"));

        TFile *dn = new TFile(dataNegName);
        effs.dataNeg.loadEff((TH2D *)dn->Get("hEffEtaPt"), (TH2D *)dn->Get("hErrlEtaPt"), (TH2D *)dn->Get("hErrhEtaPt"));

        TFile *mp = new TFile(mcPosName);
        effs.mcPos.loadEff((TH2D *)mp->Get("hEffEtaPt"), (TH2D *)mp->Get("hErrlEtaPt"), (TH2D *)mp->Get("hErrhEtaPt"));

        TFile *mn = new TFile(mcNegName);
        effs.mcNeg.loadEff((TH2D *)mn->Get("hEffEtaPt"), (TH2D *)mn->Get("hErrlEtaPt"), (TH2D *)mn->Get("hErrhEtaPt"));

        return;
    }

    void loadUncertaintyFile(effUnc &uncs, TString filename)
    {
        TFile *fSys = new TFile(filename);

        cout << "Loading ... " << filename.Data() << endl;

        uncs.fsrNeg.loadEff((TH2D *)fSys->Get("hFSRNeg"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));
        uncs.fsrPos.loadEff((TH2D *)fSys->Get("hFSRPos"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));

        uncs.mcNeg.loadEff((TH2D *)fSys->Get("hMCNeg"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));
        uncs.mcPos.loadEff((TH2D *)fSys->Get("hMCPos"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));

        uncs.bkgNeg.loadEff((TH2D *)fSys->Get("hBkgNeg"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));
        uncs.bkgPos.loadEff((TH2D *)fSys->Get("hBkgPos"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));

        uncs.tagNeg.loadEff((TH2D *)fSys->Get("hTagNeg"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));
        uncs.tagPos.loadEff((TH2D *)fSys->Get("hTagPos"), (TH2D *)fSys->Get(""), (TH2D *)fSys->Get(""));

        // cout << uncs.fsrNeg.getEff(0, 30) << endl;

        return;
    }

    double statUnc(basicEff eff, TLorentzVector *l1, int q1, TH2D *pos, TH2D *neg, double wgt, bool scale = false)
    {
        double var = 0.0;
        if (q1 > 0)
        {
            Double_t effdata = eff.dataPos.getEff(l1->Eta(), l1->Pt());
            Double_t errdata = TMath::Max(eff.dataPos.getErrLow(l1->Eta(), l1->Pt()), eff.dataPos.getErrHigh(l1->Eta(), l1->Pt()));
            // Double_t errdata = 0.5 * (eff.dataPos.getErrLow(l1->Eta(), l1->Pt()) + eff.dataPos.getErrHigh(l1->Eta(), l1->Pt()));
            Double_t effmc = eff.mcPos.getEff(l1->Eta(), l1->Pt());
            Double_t errmc = TMath::Max(eff.mcPos.getErrLow(l1->Eta(), l1->Pt()), eff.mcPos.getErrHigh(l1->Eta(), l1->Pt()));
            // Double_t errmc = 0.5* (eff.mcPos.getErrLow(l1->Eta(), l1->Pt()) + eff.mcPos.getErrHigh(l1->Eta(), l1->Pt()));
            if (scale)
            {
                // double the statistical uncertainty if applying the scale
                // because sometimes the corrections are measured charge inclusively
                // where the statistical uncertainty is underestimated by roughly 1/sqrt(2)
                errmc = errmc * sqrt(2);
                errdata = errdata * sqrt(2);
            }
            // err(a/b) = a/b * sqrt( err_a^2 / a^2 + err_b^2 / b^2)
            Double_t errSta = (effdata / effmc) * sqrt(errdata * errdata / (effdata * effdata) + errmc * errmc / (effmc * effmc));

            pos->Fill(l1->Eta(), l1->Pt(), errSta * wgt);
            var += errSta * errSta;
        }
        else
        {
            Double_t effdata = eff.dataNeg.getEff(l1->Eta(), l1->Pt());
            Double_t errdata = TMath::Max(eff.dataNeg.getErrLow(l1->Eta(), l1->Pt()), eff.dataNeg.getErrHigh(l1->Eta(), l1->Pt()));
            Double_t effmc = eff.mcNeg.getEff(l1->Eta(), l1->Pt());
            Double_t errmc = TMath::Max(eff.mcNeg.getErrLow(l1->Eta(), l1->Pt()), eff.mcNeg.getErrHigh(l1->Eta(), l1->Pt()));
            if (scale)
            {
                // similar idea as positive
                errmc *= sqrt(2);
                errdata *= sqrt(2);
            }
            Double_t errSta = (effdata / effmc) * sqrt(errdata * errdata / (effdata * effdata) + errmc * errmc / (effmc * effmc));

            neg->Fill(l1->Eta(), l1->Pt(), errSta * wgt);
            var += errSta * errSta;
        }
        return var;
    }

    double statUncHLTDilep(TLorentzVector *lep1, int q1, TLorentzVector *lep2 = nullptr, int q2 = 0)
    {
        // double effdata1 = 1.0, effdata2 = 1.0;
        // double effmc1 = 1.0, effmc2 = 1.0;
        double deff1 = 0.0, deff2 = 0.0;
        double corr1 = 1.0, corr2 = 1.0;
        // double test1 = 1.0, test2 = 1.0;
        double var = 0.0;
        TLorentzVector *l1, *l2;
        if (q1 > 0)
        {
            l1 = lep1;
            l2 = lep2;
        }
        else
        {
            l1 = lep2;
            l2 = lep1;
        }
        // sqrt((1-e2)^2*de1^2 + (1-e1)^2*de2^2))
        Double_t effdata = hlt.dataPos.getEff(l1->Eta(), l1->Pt());
        Double_t errdata = TMath::Max(hlt.dataPos.getErrLow(l1->Eta(), l1->Pt()), hlt.dataPos.getErrHigh(l1->Eta(), l1->Pt()));
        Double_t effmc = hlt.mcPos.getEff(l1->Eta(), l1->Pt());
        Double_t errmc = TMath::Max(hlt.mcPos.getErrLow(l1->Eta(), l1->Pt()), hlt.mcPos.getErrHigh(l1->Eta(), l1->Pt()));
        deff1 = sqrt(errdata * errdata / effdata / effdata + errmc * errmc / effmc / effmc);
        // double deff1d = sqrt(errdata * errdata / effdata / effdata);
        // double deff1m = sqrt(errmc * errmc / effmc / effmc);
        double deff1d = errdata;
        double deff1m = errmc;

        double corr1d = (1 - effdata);
        double corr1m = (1 - effmc);

        effdata = hlt.dataNeg.getEff(l2->Eta(), l2->Pt());
        errdata = TMath::Max(hlt.dataNeg.getErrLow(l2->Eta(), l2->Pt()), hlt.dataNeg.getErrHigh(l2->Eta(), l2->Pt()));
        effmc = hlt.mcNeg.getEff(l2->Eta(), l2->Pt());
        errmc = TMath::Max(hlt.mcNeg.getErrLow(l2->Eta(), l2->Pt()), hlt.mcNeg.getErrHigh(l2->Eta(), l2->Pt()));
        deff2 = sqrt(errdata * errdata / effdata / effdata + errmc * errmc / effmc / effmc);
        // double deff2d = sqrt(errdata * errdata / effdata / effdata);
        // double deff2m = sqrt(errmc * errmc / effmc / effmc);
        double deff2d = errdata;
        double deff2m = errmc;

        double corr2d = (1 - effdata);
        double corr2m = (1 - effmc);

        double corr = (1 - corr1d * corr2d) / (1 - corr1m * corr2m);
        // std::cout << " hahah " << corr << std::endl;

        // var = corr*sqrt(deff1*deff1+deff2*deff2);
        // var += corr1d * corr1d * deff2d * deff2d;
        // var += corr1m * corr1m * deff2m * deff2m;
        // var += corr2d * corr2d * deff1d * deff1d;
        // var += corr2m * corr2m * deff1m * deff1m;
        var += corr1d * corr1d * deff2d * deff2d;
        var += corr2d * corr2d * deff1d * deff1d;
        var += corr * corr * corr1m * corr1m * deff2m * deff2m;
        var += corr * corr * corr2m * corr2m * deff1m * deff1m;
        var = var / (1 - corr1m * corr2m) / (1 - corr1m * corr2m);

        return var;
    }

    double computeHLTSF(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        double effdata = effHLT(l1, q1, l2, q2, 1);
        double effmc = effHLT(l1, q1, l2, q2, 0);
        return effdata / effmc;
    }

    double computeStaSF(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        double effdata = effSta(l1, q1, l2, q2, 1);
        double effmc = effSta(l1, q1, l2, q2, 0);
        return effdata / effmc;
    }

    double computeSelSF(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0)
    {
        double effdata = effSel(l1, q1, l2, q2, 1);
        double effmc = effSel(l1, q1, l2, q2, 0);
        return effdata / effmc;
    }

    double effHLT(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0, bool isData = 1)
    {
        double eff = 1.0;
        if (q1 > 0)
        {
            if (isData)
                eff *= (1. - hlt.dataPos.getEff((l1->Eta()), l1->Pt()));
            else
                eff *= (1. - hlt.mcPos.getEff((l1->Eta()), l1->Pt()));
        }
        else
        {
            if (isData)
                eff *= (1. - hlt.dataNeg.getEff((l1->Eta()), l1->Pt()));
            else
                eff *= (1. - hlt.mcNeg.getEff((l1->Eta()), l1->Pt()));
        }
        if (!l2)
        {
            return 1.0 - eff;
        }
        if (q2 > 0)
        {
            if (isData)
                eff *= (1. - hlt.dataPos.getEff((l2->Eta()), l2->Pt()));
            else
                eff *= (1. - hlt.mcPos.getEff((l2->Eta()), l2->Pt()));
        }
        else
        {
            if (isData)
                eff *= (1. - hlt.dataNeg.getEff((l2->Eta()), l2->Pt()));
            else
                eff *= (1. - hlt.mcNeg.getEff((l2->Eta()), l2->Pt()));
        }
        return 1.0 - eff;
    }

    double effSta(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0, bool isData = 1)
    {
        double eff = 1.0;
        if (q1 > 0)
        {
            if (isData)
                eff *= sta.dataPos.getEff((l1->Eta()), l1->Pt());
            else
                eff *= sta.mcPos.getEff((l1->Eta()), l1->Pt());
        }
        else
        {
            if (isData)
                eff *= sta.dataNeg.getEff((l1->Eta()), l1->Pt());
            else
                eff *= sta.mcNeg.getEff((l1->Eta()), l1->Pt());
        }
        if (!l2)
        {
            return eff;
        }
        if (q2 > 0)
        {
            if (isData)
                eff *= sta.dataPos.getEff((l2->Eta()), l2->Pt());
            else
                eff *= sta.mcPos.getEff((l2->Eta()), l2->Pt());
        }
        else
        {
            if (isData)
                eff *= sta.dataNeg.getEff((l2->Eta()), l2->Pt());
            else
                eff *= sta.mcNeg.getEff((l2->Eta()), l2->Pt());
        }
        return eff;
    }

    double effSel(TLorentzVector *l1, int q1, TLorentzVector *l2 = nullptr, int q2 = 0, bool isData = 1)
    {
        double eff = 1.0;
        if (q1 > 0)
        {
            if (isData)
                eff *= sel.dataPos.getEff((l1->Eta()), l1->Pt());
            else
                eff *= sel.mcPos.getEff((l1->Eta()), l1->Pt());
        }
        else
        {
            if (isData)
                eff *= sel.dataNeg.getEff((l1->Eta()), l1->Pt());
            else
                eff *= sel.mcNeg.getEff((l1->Eta()), l1->Pt());
        }
        if (!l2)
        {
            return eff;
        }
        if (q2 > 0)
        {
            if (isData)
                eff *= sel.dataPos.getEff((l2->Eta()), l2->Pt());
            else
                eff *= sel.mcPos.getEff((l2->Eta()), l2->Pt());
        }
        else
        {
            if (isData)
                eff *= sel.dataNeg.getEff((l2->Eta()), l2->Pt());
            else
                eff *= sel.mcNeg.getEff((l2->Eta()), l2->Pt());
        }
        return eff;
    }
};

#endif
