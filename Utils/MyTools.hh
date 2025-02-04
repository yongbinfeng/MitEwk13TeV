#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <iostream>
namespace toolbox
{
    Double_t deltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2);

    Double_t deltaPhi(const Double_t phi1, const Double_t phi2);

    Int_t roundToInt(const Double_t x);

    TH1D *makeDiffHist(TH1D *hData, TH1D *hFit, const TString name);

    Int_t flavor(TClonesArray *genPartArr, Int_t vid);

    void fillGen(TClonesArray *genPartArr, Int_t vid, TLorentzVector *&vec, TLorentzVector *&lep1, TLorentzVector *&lep2, Int_t *lep1q, Int_t *lep2q, Int_t absM);

    void fillGenBorn(TClonesArray *genPartArr, Int_t vid, TLorentzVector *vec, TLorentzVector *lep1, TLorentzVector *lep2, TLorentzVector *lep3, TLorentzVector *lep4);

    float getGenLep(TClonesArray *genPartArr, TLorentzVector lep);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
    const Double_t pi = 3.14159265358979;
    Double_t dphi = fabs(phi1 - phi2);
    while (dphi > pi)
        dphi = fabs(dphi - 2.0 * pi);

    Double_t deta = eta1 - eta2;

    return sqrt(dphi * dphi + deta * deta);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaPhi(Double_t phi1, Double_t phi2)
{
    // Compute dPhi between two given angles. Results is in [0,pi].
    const Double_t pi = 3.14159265358979;
    Double_t dphi = fabs(phi1 - phi2);
    while (dphi > pi)
        dphi = fabs(dphi - 2.0 * pi);

    return dphi;
}

//------------------------------------------------------------------------------------------------------------------------
Int_t toolbox::roundToInt(Double_t x)
{
    if (x > 0)
        return ((x - floor(x)) < (ceil(x) - x)) ? (Int_t)floor(x) : (Int_t)ceil(x);
    else
        return ((x - floor(x)) < (ceil(x) - x)) ? (Int_t)ceil(x) : (Int_t)floor(x);
}

TH1D *toolbox::makeDiffHist(TH1D *hData, TH1D *hFit, const TString name)
{
    TH1D *hDiff = (TH1D *)hData->Clone("hDiff");
    hDiff->SetName(name);
    for (Int_t ibin = 1; ibin <= hData->GetNbinsX(); ibin++)
    {

        Double_t diff0 = (hData->GetBinContent(ibin) - hFit->GetBinContent(ibin));
        Double_t diff = 0;
        Double_t err = 0;
        if (hData->GetBinContent(ibin) != 0)
        {
            diff = diff0 / hData->GetBinContent(ibin);
            err = (hFit->GetBinContent(ibin) / hData->GetBinContent(ibin)) * sqrt((1.0 / hFit->GetBinContent(ibin)) + (1.0 / hData->GetBinContent(ibin)));
        }
        hDiff->SetBinContent(ibin, diff);
        hDiff->SetBinError(ibin, err);
    }

    hDiff->GetYaxis()->SetTitleOffset(0.42);
    hDiff->GetYaxis()->SetTitleSize(0.13);
    hDiff->GetYaxis()->SetLabelSize(0.10);
    hDiff->GetYaxis()->SetNdivisions(104);
    hDiff->GetYaxis()->CenterTitle();
    hDiff->GetXaxis()->SetTitleOffset(1.2);
    hDiff->GetXaxis()->SetTitleSize(0.13);
    hDiff->GetXaxis()->SetLabelSize(0.12);
    hDiff->GetXaxis()->CenterTitle();

    return hDiff;
}

//------------------------------------------------------------------------------------------------------------------------
Int_t toolbox::flavor(TClonesArray *genPartArr, Int_t vid)
{

    for (Int_t i = 0; i < genPartArr->GetEntries(); i++)
    {
        const baconhep::TGenParticle *genloop = (baconhep::TGenParticle *)((*genPartArr)[i]);
        Int_t pdgId = fabs(genloop->pdgId);
        Int_t parentPdgId = fabs(dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent > -1 ? genloop->parent : 0))->pdgId);
        if ((pdgId == 11 || pdgId == 13 || pdgId == 15) && (parentPdgId == fabs(vid) || genloop->status == 23))
        {
            return genloop->pdgId;
        }
    }
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------
void toolbox::fillGen(TClonesArray *genPartArr, Int_t vid, TLorentzVector *&vec, TLorentzVector *&lep1, TLorentzVector *&lep2, Int_t *lep1q, Int_t *lep2q, Int_t absM)
{
    Int_t iv = -1, iv1 = -1, iv2 = -1;
    TLorentzVector *lepPos = 0, *lepNeg = 0;
    TLorentzVector *preLepPos = 0, *preLepNeg = 0;
    for (Int_t i = 0; i < genPartArr->GetEntries(); i++)
    {
        const baconhep::TGenParticle *genloop = (baconhep::TGenParticle *)((*genPartArr)[i]);
        // if (fabs(genloop->pdgId) != 22)
        //     std::cout << i << " " << genloop->pdgId << " " << genloop->status << " " << genloop->parent << " " <<dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>0 ? genloop->parent : 1))->pdgId << " " << genloop->pt << " " << genloop->mass << std::endl;
        //     else if (fabs(genloop->pdgId) == 22 && fabs(dynamic_cast<baconhep::TGenParticle*>(genPartArr->At(genloop->parent > 0 ? genloop->parent : 1))->pdgId) == 11)
        //             std::cout
        //         << i << " " << genloop->pdgId << " " << genloop->status << " " << genloop->parent << " " << dynamic_cast<baconhep::TGenParticle*>(genPartArr->At(genloop->parent > 0 ? genloop->parent : 1))->pdgId << " " << genloop->pt << " " << genloop->mass << std::endl;

        if (genloop->status == 44 && (fabs(genloop->pdgId) == 15 || fabs(genloop->pdgId) == 13 || fabs(genloop->pdgId) == 11 || fabs(genloop->pdgId) == 16 || fabs(genloop->pdgId) == 14 || fabs(genloop->pdgId) == 12))
        {
            if (genloop->pdgId < 0 && lepPos == 0)
            {
                lepPos = new TLorentzVector(0, 0, 0, 0);
                lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                preLepPos = new TLorentzVector(0, 0, 0, 0);
                preLepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                iv1 = i;
            }
            else if (genloop->pdgId > 0 && lepNeg == 0)
            {
                lepNeg = new TLorentzVector(0, 0, 0, 0);
                lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                preLepNeg = new TLorentzVector(0, 0, 0, 0);
                preLepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                iv2 = i;
            }
        }
        else if ((absM == 0 && genloop->pdgId == vid && (genloop->status == 3 || genloop->status == 22)) || (absM == 1 && fabs(genloop->pdgId) == fabs(vid) && (genloop->status == 3 || genloop->status == 62)))
        {
            vec = new TLorentzVector(0, 0, 0, 0);
            vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv = i;
        }
        else if (iv != -1 && genloop->parent == iv)
        {
            if ((absM == 0 && genloop->pdgId == vid) || (absM == 1 && fabs(genloop->pdgId) == fabs(vid)))
            {
                vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                iv = i;
            }
            else if (fabs(genloop->pdgId) == 15 || fabs(genloop->pdgId) == 13 || fabs(genloop->pdgId) == 11 || fabs(genloop->pdgId) == 16 || fabs(genloop->pdgId) == 14 || fabs(genloop->pdgId) == 12)
            {
                if (genloop->pdgId < 0 && lepPos == 0)
                {
                    lepPos = new TLorentzVector(0, 0, 0, 0);
                    lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                    preLepPos = new TLorentzVector(0, 0, 0, 0);
                    preLepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                    iv1 = i;
                }
                else if (genloop->pdgId > 0 && lepNeg == 0)
                {
                    lepNeg = new TLorentzVector(0, 0, 0, 0);
                    lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                    preLepNeg = new TLorentzVector(0, 0, 0, 0);
                    preLepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                    iv2 = i;
                }
            }
        }
        else if (iv1 != -1 && genloop->parent == iv1)
        {
            lepPos->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv1 = i;
        }
        else if (iv2 != -1 && genloop->parent == iv2)
        {
            lepNeg->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv2 = i;
        }
    }

    if (iv == -1 && preLepNeg && preLepPos)
    {
        TLorentzVector temp = *preLepNeg + *preLepPos;
        vec->SetPtEtaPhiM(temp.Pt(), temp.Eta(), temp.Phi(), temp.M());
        // vec->SetPtEtaPhiM(0,0,0,temp.M());
        // std::cout << "No Boson " << preLepNeg->Pt() << " " << preLepPos->Pt() << " " << temp.M() << " " << temp.Pt() << " " << temp.Eta() << " " << temp.Phi() << " " <<  vec->M() << std::endl;
    }

    if (lepPos && lepNeg)
    {
        if (lepPos->Pt() > lepNeg->Pt())
        {
            lep1->SetPtEtaPhiM(lepPos->Pt(), lepPos->Eta(), lepPos->Phi(), lepPos->M());
            lep2->SetPtEtaPhiM(lepNeg->Pt(), lepNeg->Eta(), lepNeg->Phi(), lepNeg->M());
            *lep1q = 1;
            *lep2q = -1;
        }
        else
        {
            lep1->SetPtEtaPhiM(lepNeg->Pt(), lepNeg->Eta(), lepNeg->Phi(), lepNeg->M());
            lep2->SetPtEtaPhiM(lepPos->Pt(), lepPos->Eta(), lepPos->Phi(), lepPos->M());
            *lep1q = -1;
            *lep2q = 1;
        }
    }
    else if (lepPos)
    {
        lep1->SetPtEtaPhiM(lepPos->Pt(), lepPos->Eta(), lepPos->Phi(), lepPos->M());
        *lep1q = 1;
        *lep2q = 0;
    }
    else if (lepNeg)
    {
        lep1->SetPtEtaPhiM(lepNeg->Pt(), lepNeg->Eta(), lepNeg->Phi(), lepNeg->M());
        *lep1q = -1;
        *lep2q = 0;
    }
    // std::cout << "Vector boson " << vec->Pt() << " " << vec->Eta() << std::endl;
    delete preLepNeg;
    delete preLepPos;
    delete lepNeg;
    delete lepPos;

    preLepNeg = 0;
    preLepPos = 0;
    lepNeg = 0;
    lepPos = 0;
}

//------------------------------------------------------------------------------------------------------------------------
void toolbox::fillGenBorn(TClonesArray *genPartArr, Int_t vid, TLorentzVector *vec, TLorentzVector *lep1, TLorentzVector *lep2, TLorentzVector *lep3, TLorentzVector *lep4)
{
    // lep1, lep3 are the paricles, lep2 and lep4 are the anti-particles
    // lep1,lep2->pre-fsr. lep3,lep4->post-fsr
    Int_t iv = -1, iv1 = -1, iv2 = -1;
    Int_t pdgId1 = -1, pdgId2 = -1;
    for (Int_t i = 0; i < genPartArr->GetEntries(); i++)
    {
        const baconhep::TGenParticle *genloop = (baconhep::TGenParticle *)((*genPartArr)[i]);
        if (fabs(genloop->pdgId) == 22)
            continue;
        // if(fabs(genloop->pdgId)!=22)
        //  std::cout << i << " " << genloop->pdgId << " " << genloop->status << " " << genloop->parent << " " <<dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>0 ? genloop->parent : 1))->pdgId << " " << genloop->pt << " " << genloop->mass << std::endl;
        //  else if(fabs(genloop->pdgId)==22 && fabs(dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>0 ? genloop->parent : 1))->pdgId)==11)
        //  std::cout << i << " " << genloop->pdgId << " " << genloop->status << " " << genloop->parent << " " <<dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genloop->parent>0 ? genloop->parent : 1))->pdgId << " " << genloop->pt << " " << genloop->mass << std::endl;

        if (fabs(genloop->pdgId) == vid)
        {
            vec->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv = i;
        }
        else if ((genloop->pdgId == 11 || genloop->pdgId == 13 || genloop->pdgId == 12 || genloop->pdgId == 14 || genloop->pdgId == 15 || genloop->pdgId == 16) && genloop->parent == iv)
        {
            lep1->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            lep3->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv1 = i;
            pdgId1 = genloop->pdgId;
        }
        else if ((genloop->pdgId == -11 || genloop->pdgId == -13 || genloop->pdgId == -12 || genloop->pdgId == -14 || genloop->pdgId == -15 || genloop->pdgId == -16) && genloop->parent == iv)
        {
            lep2->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            lep4->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv2 = i;
            pdgId2 = genloop->pdgId;
        }
        else if ((genloop->pdgId == 11 || genloop->pdgId == 13 || genloop->pdgId == 12 || genloop->pdgId == 14 || genloop->pdgId == 15 || genloop->pdgId == 16) && iv == -1 && (genloop->status == 44 || genloop->status == 23)) // status < 50 for pythia8 means particles produced by initial state showers
        {
            lep1->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            lep3->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv1 = i;
            pdgId1 = genloop->pdgId;
        }
        else if ((genloop->pdgId == -11 || genloop->pdgId == -13 || genloop->pdgId == -12 || genloop->pdgId == -14 || genloop->pdgId == -15 || genloop->pdgId == -16) && iv == -1 && (genloop->status == 44 || genloop->status == 23)) // status < 50 for pythia8 means particles produced by initial state showers
        {
            lep2->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            lep4->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv2 = i;
            pdgId2 = genloop->pdgId;
        }
        else if (iv1 != -1 && genloop->parent == iv1 && genloop->pdgId == pdgId1)
        {
            // trace the lepton status change
            lep3->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv1 = i;
        }
        else if (iv2 != -1 && genloop->parent == iv2 && genloop->pdgId == pdgId2)
        {
            // trace the lepton status change
            lep4->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            iv2 = i;
        }
        else if (iv1 != -1 && genloop->parent == iv1 && pdgId1 == 15)
        {
            // for tau decays
            // 15 -> 11/13 + -12/-14 + 16
            if (fabs(pdgId2) == 16 && (fabs(genloop->pdgId) == 12 || fabs(genloop->pdgId) == 14 || fabs(genloop->pdgId) == 16))
            {
                // add the other neutrinos to lep2
                TLorentzVector nu;
                nu.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                (*lep2) = (*lep2) + nu;
                (*lep4) = (*lep4) + nu;
            }
            else if (fabs(genloop->pdgId) == 11 || fabs(genloop->pdgId) == 13)
            {
                // lepton decay product of the tau lepton: either muon or electron
                // this lepton could further FSR?
                lep3->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            }
        }
        else if (iv2 != -1 && genloop->parent == iv2 && pdgId2 == -15)
        {
            // for tau decays
            // -15 -> -11/-13 + 12/14 - 16
            if (fabs(pdgId1) == 16 && (fabs(genloop->pdgId) == 12 || fabs(genloop->pdgId) == 14 || fabs(genloop->pdgId) == 16))
            {
                TLorentzVector nu;
                nu.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
                (*lep1) = (*lep1) + nu;
                (*lep3) = (*lep3) + nu;
            }
            else if (fabs(genloop->pdgId) == 11 || fabs(genloop->pdgId) == 13)
            {
                // lepton decay product of the tau lepton: either muon or electron
                // this lepton could further FSR?
                lep4->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            }
        }
    }
    if (iv == -1)
        *vec = (*lep1) + (*lep2);
}

// get gen lepton invo
float toolbox::getGenLep(TClonesArray *genPartArr, TLorentzVector lep)
{
    // Int_t iv=-1, iv1=-1, iv2=-1;
    TLorentzVector *lepcand = 0;
    TLorentzVector *bestlep = 0;
    double genPt = 0;
    double bestCandDist = 9999.;
    // std::cout << "hello" << std::endl;
    // // TLorentzVector *preLepPos=0, *preLepNeg=0;
    // std::cout << "n entried in gen part array " << genPartArr->GetEntries() << std::endl;
    for (Int_t i = 0; i < genPartArr->GetEntries(); i++)
    {
        const baconhep::TGenParticle *genloop = (baconhep::TGenParticle *)((*genPartArr)[i]);
        // std::cout << "blah " << std::endl;
        // if(fabs(genloop->pdgId)==22) continue;
        // if(fabs(genloop->pdgId)!=13&&fabs(genloop->pdgId)!=11&&fabs(genloop->pdgId)!=15&&fabs(genloop->pdgId)!=12&&fabs(genloop->pdgId)!=14&&fabs(genloop->pdgId)!=16&&fabs(genloop->pdgId)!=24&&fabs(genloop->pdgId)!=23) continue;
        if (fabs(genloop->pdgId) != 13)
            continue;
        // std::cout << "okayy" << std::endl;
        // std::cout << i << " " << genloop->pdgId << " " << genloop->status << " " << genloop->parent << " " << genloop->pt << " " << genloop->mass << std::endl;

        lepcand = new TLorentzVector(0, 0, 0, 0);
        lepcand->SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
        // std::cout << "asdlfkjlkj" << std::endl;
        double lepCandDist = deltaR(lepcand->Eta(), lepcand->Phi(), lep.Eta(), lep.Phi());
        // if(deltaR(lepcand->Eta(), lepcand->Phi(), lep.Eta(), lep.Phi())>0.5) continue;
        if (genloop->status != 1)
            continue;
        // std::cout << "found candiate  " << i << "  parent_i " << genloop->parent << " status " << genloop->status << " type " << genloop->pdgId << std::endl;
        // std::cout << "gen pt " << genloop->pt << "  lep pt " << lep.Pt() << std::endl;
        // if(genloop->parent<0) continue;
        // const baconhep::TGenParticle* genparent = (baconhep::TGenParticle*) ((*genPartArr)[genloop->parent]);
        // std::cout << "gen parent info " << genparent->pdgId << " gen parent pT " << genparent->pt << std::endl;
        // if(abs(genloop->pt-lep.Pt()) < abs(genPt - lep.Pt()))
        // std::cout << "best " << bestCandDist << "  lepCand " << lepCandDist << std::endl;
        if (lepCandDist < bestCandDist)
        {
            genPt = genloop->pt;
            bestlep = lepcand;
            bestCandDist = lepCandDist;
        }
        // std::cout << "new gen pt " << genPt << std::endl;
        /// check if it's status 33 or 23
        // if(genloop->status%10==3&&)
    }
    delete lepcand;

    return genPt;

    // //std::cout << "Vector boson " << vec->Pt() << " " << vec->Eta() << std::endl;
    // delete preLepNeg;
    // delete preLepPos;
    // delete lepNeg;
    // delete lepPos;

    // preLepNeg=0; preLepPos=0;
    // lepNeg=0; lepPos=0;
}
#endif
