#define TUPLE_NEW

#ifndef DERIVEDTUPLEVARIABLES_H
#define DERIVEDTUPLEVARIABLES_H

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TopTagger/Tools/cpp/TaggerUtility.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"
#include "ScaleFactors.h"
#include "ScaleFactorsttBar.h"

#include "TopTagger.h"
#include "TTModule.h"
#include "TopTaggerUtilities.h"
#include "TopTaggerResults.h"

#include "TopTagger/TopTagger/include/TopObject.h"
#include "TopTagger/CfgParser/include/TTException.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include "TVector2.h"

#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <random>
 
namespace plotterFunctions
{
    class GenerateWeight
    {
    private:
        TH1* muEff;
        TH1* muEffReco;
        TH2* muEffIso;
        TH2* elecEffReco;
        TH2* elecEffIso;
        TH2* muAcc;
        TH1* hZEff;
        TH1* hZAcc;
        TH1* hZAccElec;
        TH1* hZAccLepPt;
        TH1* hZAccLepPtElec;

        void generateWeight(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& jetsLVec         = tr.getVec<TLorentzVector>("jetsLVecLepCleaned");
            const std::vector<TLorentzVector>& cutMuVec         = tr.getVec<TLorentzVector>("cutMuVec");
            const std::vector<TLorentzVector>& cutElecVec       = tr.getVec<TLorentzVector>("cutElecVec");
            const std::vector<float>& cutMuActivity            = tr.getVec<float>("cutMuActivity");
            const std::vector<float>& cutElecActivity          = tr.getVec<float>("cutElecActivity");
            const std::vector<TLorentzVector*>& genMu           = tr.getVec<TLorentzVector*>("genMu");
            const std::vector<TLorentzVector*>& genMuInAcc      = tr.getVec<TLorentzVector*>("genMuInAcc");
            const std::vector<TLorentzVector*>& genMatchMuInAcc = tr.getVec<TLorentzVector*>("genMatchMuInAcc");

            const int& pdgIdZDec      = tr.getVar<int>("pdgIdZDec");
            const bool& passMuZinvSel = tr.getVar<bool>("passMuZinvSel");
            const bool& passElecZinvSel = tr.getVar<bool>("passElecZinvSel");
            const float& ht          = tr.getVar<float>("ht");
            const float& bestRecoZPt = tr.getVar<float>("bestRecoZPt");
            const float& genZPt      = tr.getVar<float>("genZPt");

            const int& nJets     =  tr.getVar<int>("nJets");
	    const float& stored_weight = tr.getVar<float>("stored_weight");

            // Calculate PU weight

            // Calculate Z-eff weight

            const float zMassMin = 71.0;
            const float zMass    = 91.0;
            const float zMassMax = 111.0;

            float zMassCurrent = 1.0e300, zEff = 1.0e100, zAcc = 1.0e100;
            for(int i = 0; i < cutMuVec.size(); ++i)
            {
                if(cutMuVec[i].Pt() < 10) continue;
                for(int j = 0; j < i && j < cutMuVec.size(); ++j)
                {
                    if(cutMuVec[j].Pt() < 10) continue;
                    float zm = (cutMuVec[i] + cutMuVec[j]).M();
                    if(fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    {
                        zMassCurrent = zm;
                        float mu1pt = cutMuVec[i].Pt();
                        float mu2pt = cutMuVec[j].Pt();
                        float mu1Act = cutMuActivity[i];
                        float mu2Act = cutMuActivity[j];

                        //set to not overflow histograms
                        if(mu1pt >= 2000.0) mu1pt = 1999.9;
                        if(mu2pt >= 2000.0) mu2pt = 1999.9;

                        //Get mu efficiencies
                        float muEff1 = 0.0, muEff2 = 0.0;

                        if(muEff && muEffReco && muEffIso)
                        {
                            //Fit to reco eff (eff = p0 + p1*pt + p2*pt^2
                            //PHYS14
                            //const float fitStart = 200.0; // extended to 1400 GeV
                            //const float p0 =     0.955847; // +/- 0.461944
                            //const float p1 = -2.24431e-05; // +/- 0.00128305
                            //const float p2 = -5.68907e-08; // +/- 7.85913e-07
                            //Sprint15
                            //const float fitStart = 200.0; // extended to 1400 GeV
                            //const float p0 =  9.83467e-01; // +/- 1.54469e+00
                            //const float p1 = -7.81897e-06; // +/- 4.16487e-03
                            //const float p2 = -1.22092e-08; // +/- 2.18556e-06
                            //Spring15 low stats
                            //const float fitStart = 200.0;  // extended to 2000 GeV
                            //const float p0 =     0.979238; //   +/-   1.17313
                            //const float p1 = -6.47338e-06; //   +/-   0.00342715
                            //const float p2 = -1.16258e-08; //   +/-   1.87837e-06
                            //Spring15 extended samples
                            const float fitStart = 200.0;  // extended to 2000 GeV
                            const float p0 =      0.97431; //   +/-   1.17119     
                            const float p1 =  2.87484e-05; //   +/-   0.00341371  
                            const float p2 = -5.37058e-08; //   +/-   1.86309e-06 

                            float muRecoEff = 0.0;

                            int recoPtBin = muEffReco->GetXaxis()->FindBin(mu1pt);
                            if(mu1pt > fitStart) muRecoEff = p0 + p1*mu1pt + p2*mu1pt*mu1pt;
                            else                 muRecoEff = muEffReco->GetBinContent(recoPtBin);
                            int isoPtBin = muEffIso->GetXaxis()->FindBin(mu1pt);
                            if(isoPtBin >= muEffIso->GetNbinsX()) isoPtBin = muEffIso->GetNbinsX();
                            int isoActBin = muEffIso->GetYaxis()->FindBin(mu1Act);
                            if(isoActBin >= muEffIso->GetNbinsY()) isoActBin = muEffIso->GetNbinsY();
                            muEff1 = muRecoEff * muEffIso->GetBinContent(isoPtBin, isoActBin);

                            muRecoEff = 0.0;
                            recoPtBin = muEffReco->GetXaxis()->FindBin(mu2pt);
                            if(mu2pt > fitStart) muRecoEff = p0 + p1*mu2pt + p2*mu2pt*mu2pt;
                            else                 muRecoEff = muEffReco->GetBinContent(recoPtBin);
                            isoPtBin = muEffIso->GetXaxis()->FindBin(mu2pt);
                            if(isoPtBin >= muEffIso->GetNbinsX()) isoPtBin = muEffIso->GetNbinsX();
                            isoActBin = muEffIso->GetYaxis()->FindBin(mu2Act);
                            if(isoActBin >= muEffIso->GetNbinsY()) isoActBin = muEffIso->GetNbinsY();
                            muEff2 = muRecoEff * muEffIso->GetBinContent(isoPtBin, isoActBin);
                        }


                        if((mu1pt > 20 && muEff1 < 1.0e-5) || (mu2pt > 20 && muEff2 < 1.0e-5))
                        {
                            std::cout << "SMALL muEff!!! muEff1: " << muEff1 << "\tmuEff2: " << muEff2 << "\t" << mu1pt << "\t" << mu2pt << std::endl;
                            zEff = 1.0e-10;
                        }
                        else zEff = muEff1 * muEff2;
                    }
                }
            }

            zMassCurrent = 1.0e300;
            float zEffElec = 1.0e100, zAccElec = 1.0e100;
            for(int i = 0; i < cutElecVec.size(); ++i)
            {
                if(cutElecVec[i].Pt() < 10) continue;
                for(int j = 0; j < i && j < cutElecVec.size(); ++j)
                {
                    if(cutElecVec[j].Pt() < 10) continue;
                    float zm = (cutElecVec[i] + cutElecVec[j]).M();
                    if(fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    {
                        zMassCurrent = zm;
                        float elec1pt = cutElecVec[i].Pt();
                        float elec2pt = cutElecVec[j].Pt();

                        float elec1Act = cutElecActivity[i];
                        float elec2Act = cutElecActivity[j];

                        //set to not overflow histograms
                        if(elec1pt >= 2000.0) elec1pt = 1999.9;
                        if(elec2pt >= 2000.0) elec2pt = 1999.9;

                        //Get elec efficiencies
                        float elecEff1 = 0.0, elecEff2 = 0.0;
                        if(elecEffReco && elecEffIso)
                        {
                            //Fit to iso eff (eff = p0 + p1*pt + p2*pt^2
                            const float fitStart =200.0;  // extended to 2000 GeV
                            const float p0 =     0.989411; //   +/-   1.18467
                            const float p1 =  3.66321e-06; //   +/-   0.00346729
                            const float p2 = -5.68292e-09; //   +/-   1.90334e-06

                            float elecIsoEff = 0.0;

                            //int isoPtBin = elecEffIso->GetXaxis()->FindBin(elec1pt);
                            //if(elec1pt > fitStart) elecIsoEff = p0 + p1*elec1pt + p2*elec1pt*elec1pt;
                            //else                   elecIsoEff = elecEffIso->GetBinContent(isoPtBin);

                            int recoPtBin = elecEffReco->GetXaxis()->FindBin(elec1pt);
                            if(recoPtBin >= elecEffReco->GetNbinsX()) recoPtBin = elecEffReco->GetNbinsX();

                            int recoActBin = elecEffReco->GetYaxis()->FindBin(elec1Act);
                            if(recoActBin >= elecEffReco->GetNbinsY()) recoActBin = elecEffReco->GetNbinsY();

                            elecIsoEff = elecEffIso->GetBinContent(recoPtBin, recoActBin);
                            elecEff1 = elecIsoEff * elecEffReco->GetBinContent(recoPtBin, recoActBin);
                            //std::cout << elec1pt << "\t" << recoPtBin << "\t" <<  elec1Act << "\t" << recoActBin << "\t" << elecEff1 << "\t" << elecIsoEff << "\t" << elecEffReco->GetBinContent(recoPtBin, recoActBin) << std::endl;

                            elecIsoEff = 0.0;
                            //isoPtBin = elecEffIso->GetXaxis()->FindBin(elec2pt);
                            //if(elec2pt > fitStart) elecIsoEff = p0 + p1*elec2pt + p2*elec2pt*elec2pt;
                            //else                   elecIsoEff = elecEffIso->GetBinContent(isoPtBin);

                            recoPtBin = elecEffReco->GetXaxis()->FindBin(elec2pt);
                            if(recoPtBin >= elecEffReco->GetNbinsX()) recoPtBin = elecEffReco->GetNbinsX();

                            recoActBin = elecEffReco->GetYaxis()->FindBin(elec2Act);
                            if(recoActBin >= elecEffReco->GetNbinsY()) recoActBin = elecEffReco->GetNbinsY();

                            elecIsoEff = elecEffIso->GetBinContent(recoPtBin, recoActBin);
                            elecEff2 = elecIsoEff * elecEffReco->GetBinContent(recoPtBin, recoActBin);
                            //std::cout << elec2pt << "\t" << recoPtBin << "\t" << elec2Act << "\t" << recoActBin << "\t" << elecEff2 << "\t" << elecIsoEff << "\t" << elecEffReco->GetBinContent(recoPtBin, recoActBin) << std::endl;
                        }


                        if((elec1pt > 33 && elec2pt) && (elecEff1 < 1.0e-5 || elecEff2 < 1.0e-5))
                        {
                            std::cout << "SMALL elecEff!!! elecEff1: " << elecEff1 << "\telecEff2: " << elecEff2 << "\t" << elec1pt << "\t" << elec2pt << std::endl;
                            zEffElec = 1.0e-10;
                        }
                        else zEffElec = elecEff1 * elecEff2;
                    }
                }
            }

            float genCleanHt = ht;
            for(auto& tlvp : genMuInAcc) if(tlvp->Pt() > 50) genCleanHt -= tlvp->Pt();

            const float MHT_jetPtMin = 30.0;
            TLorentzVector MHT;
            float mu1dRMin = 99.9, mu2dRMin = 99.9;
            for(auto& jet : jetsLVec)
            {
                float mu1dR = 999.9, mu2dR = 999.9;
                if(cutMuVec.size() >= 1) mu1dR = ROOT::Math::VectorUtil::DeltaR(jet, cutMuVec[0]);
                if(cutMuVec.size() >= 2) mu2dR = ROOT::Math::VectorUtil::DeltaR(jet, cutMuVec[1]);
                mu1dRMin = std::min(mu1dRMin, mu1dR);
                mu2dRMin = std::min(mu2dRMin, mu2dR);

                if(jet.Pt() > MHT_jetPtMin) MHT += jet;
            }

            float mudR = 99.9, genMudR = 99.9, genMudPhi = 99.9, genMudEta = 99.9;

            if(cutMuVec.size() >= 2) mudR    = ROOT::Math::VectorUtil::DeltaR(cutMuVec[0], cutMuVec[1]);
            if(genMu.size() > 2) std::cout << "MORE THAN 2 GEN MUONS: " << genMu.size() << std::endl;
            if(genMu.size() >= 2)
            {
                genMudR = ROOT::Math::VectorUtil::DeltaR(*(genMu[0]), *(genMu[1]));
                genMudPhi = ROOT::Math::VectorUtil::DeltaPhi(*(genMu[0]), *(genMu[1]));
                genMudEta = genMu[0]->Eta() - genMu[1]->Eta();
            }
            //std::cout<<"cutMuVec at size of two "<<cutMuVec.size()<<std::endl;
            // muon Z pt acceptance corrections
            //functional form [2] - exp([0] + [1]*x)
            //PHYS14
            //float acc_p0 = -2.91374e-01;
            //float acc_p1 = -4.42884e-03;
            //float acc_p2 =  9.51190e-01;
            //Spring15 low stats
            //float acc_p0 = -2.64921e-01;
            //float acc_p1 = -4.65305e-03;
            //float acc_p2 =  9.50493e-01;
            //Spring15 extended sample
            const float acc_p0 = -2.60663e-01;
            const float acc_p1 = -4.60623e-03;
            const float acc_p2 = 9.49434e-01;

            if(hZAcc && hZAccLepPt)
            {
                //Eta portion of acceptance
                if(bestRecoZPt < 150) zAcc = hZAcc->GetBinContent(hZAcc->GetXaxis()->FindBin(bestRecoZPt));
                else                  zAcc = acc_p2 - exp(acc_p0 + acc_p1 * bestRecoZPt);

                //pt portion of acceptance
                zAcc *= hZAccLepPt->GetBinContent(hZAccLepPt->GetXaxis()->FindBin(bestRecoZPt));
            }

            if(passMuZinvSel && zAcc < 0.05)
            {
                std::cout << "WARNING: Muon Z acceptance < 0.05, forcing weight to zero! Acc_Z: " << zAcc << std::endl;
                zAcc = 1.0e101;
            }

            if(passMuZinvSel && zEff < 0.04)
            {
                std::cout << "WARNING: Muon Z efficiency < 0.05, forcing weight to zero! Eff_Z: " << zEff << "\tZ(Pt): " << bestRecoZPt <<  std::endl;
                zEff = 1.0e101;
            }

            // electron Z pt acceptance corrections
            //functional form [2] - exp([0] + [1]*x)
            //Spring15
            //float elecAcc_p0 = -3.18955e-02;
            //float elecAcc_p1 = -5.16522e-03;
            //float elecAcc_p2 = 8.92280e-01;
            float elecAcc_p0 = -2.50679e-01;
            float elecAcc_p1 = -5.34976e-03;
            float elecAcc_p2 = 9.33562e-01;

            if(hZAccElec && hZAccLepPtElec)
            {
                //Eta portion of the acceptance
                if(bestRecoZPt < 100) zAccElec = hZAccElec->GetBinContent(hZAccElec->GetXaxis()->FindBin(bestRecoZPt));
                else                  zAccElec = elecAcc_p2 - exp(elecAcc_p0 + elecAcc_p1 * bestRecoZPt);

                //Pt portion of the acceptance
                zAccElec = hZAccLepPtElec->GetBinContent(hZAccLepPtElec->GetXaxis()->FindBin(bestRecoZPt));
            }

            if(passElecZinvSel && zAccElec < 0.05)
            {
                std::cout << "WARNING: Elec Z acceptance < 0.05, forcing weight to zero! Acc_Z: " << zAccElec << std::endl;
                zAccElec = 1.0e101;
            }

            if(passElecZinvSel && zEffElec < 0.04)
            {
                std::cout << "WARNING: Elec Z efficiency < 0.05, forcing weight to zero! Eff_Z: " << zEffElec << "\tZ(Pt): " << bestRecoZPt <<  std::endl;
                zEffElec = 1.0e101;
            }

	    // Process the generator weight
	    float genWeight = 1.;
	    // Never apply this weight for data! In the old ntuple version <=3 this is "-1", in the newer ones it is "0"
	    if(stored_weight < 0)
	      genWeight = -1.;

            tr.registerDerivedVar("mu1dRMin", mu1dRMin);
            tr.registerDerivedVar("mu2dRMin", mu2dRMin);
            tr.registerDerivedVar("mudR", mudR);
            tr.registerDerivedVar("genMudR", genMudR);
            tr.registerDerivedVar("genMudPhi", genMudPhi);
            tr.registerDerivedVar("genMudEta", genMudEta);

            tr.registerDerivedVar("zEff", zEff);
            tr.registerDerivedVar("zEffWgt", 1.0/zEff);
            tr.registerDerivedVar("zAcc", zAcc);
            tr.registerDerivedVar("zAccWgt", 1.0/zAcc);

            tr.registerDerivedVar("zEffElec", zEffElec);
            tr.registerDerivedVar("zEffWgtElec", 1.0/zEffElec);
            tr.registerDerivedVar("zAccElec", zAccElec);
            tr.registerDerivedVar("zAccWgtElec", 1.0/zAccElec);

            tr.registerDerivedVar("cleanMHt", MHT.Pt());
            tr.registerDerivedVar("cleanMHtPhi", MHT.Phi());

	    tr.registerDerivedVar("genWeight", genWeight);
        }

    public:
        GenerateWeight()
        {
            muEff          = nullptr;
            muEffReco      = nullptr;
            muEffIso       = nullptr;
            muAcc          = nullptr;
            hZEff          = nullptr;
            hZAcc          = nullptr;
            elecEffReco    = nullptr;
            elecEffIso     = nullptr;
            hZAccElec      = nullptr;
            hZAccLepPt     = nullptr;
            hZAccLepPtElec = nullptr;

            TH1::AddDirectory(false);
            TFile *f = new TFile("lepEffHists.root");
            if(f)
            {
                muEff          = static_cast<TH1*>(f->Get("hMuEffPt_ratio"));
                muEffReco      = static_cast<TH1*>(f->Get("hMuEffPtReco_ratio"));
                muEffIso       = static_cast<TH2*>(f->Get("hMuEffPtActIso_ratio"));
                muAcc          = static_cast<TH2*>(f->Get("hMuAcc_ratio"));
                hZEff          = static_cast<TH1*>(f->Get("hZEffPt_ratio"));
                hZAcc          = static_cast<TH1*>(f->Get("hZAccPtSmear_ratio"));
                elecEffReco    = static_cast<TH2*>(f->Get("hElecEffPtActReco_ratio"));
                elecEffIso     = static_cast<TH2*>(f->Get("hElecEffPtActIso_ratio"));
                hZAccElec      = static_cast<TH1*>(f->Get("hZElecAccPtSmear_ratio"));
                hZAccLepPt     = static_cast<TH1*>(f->Get("hZAccPtMuPtSmear_ratio"));
                hZAccLepPtElec = static_cast<TH1*>(f->Get("hZElecAccPtPtSmear_ratio"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open: lepEffHists.root" << std::endl;
            }
        }

        ~GenerateWeight()
        {
            //if(muEff)          delete muEff;
            //if(muEffReco)      delete muEffReco;
            //if(muEffIso)       delete muEffIso;
            //if(muAcc)          delete muAcc;
            //if(hZEff)          delete hZEff;
            //if(hZAcc)          delete hZAcc;
            //if(elecEffReco)    delete elecEffReco;
            //if(elecEffIso)     delete elecEffIso;
            //if(hZAccElec)      delete hZAccElec;
            //if(hZAccLepPt)     delete hZAccLepPt;
            //if(hZAccLepPtElec) delete hZAccLepPtElec;
        }

        void operator()(NTupleReader& tr)
        {
            generateWeight(tr);
        }
    };

    class NJetWeight
    {
    private:
        TH1* njWTTbar_0b;
        TH1* njWDYZ_0b;
        TH1* njWTTbar_g1b;
        TH1* njWDYZ_g1b;

        TH1* MCfake1b;
        TH1* MCfake2b;
        TH1* MCfake3b;

        void generateWeight(NTupleReader& tr)
        {
            const int& cntNJetsPt30Eta24Zinv = tr.getVar<int>("cntNJetsPt30Eta24Zinv");
            const int& cntCSVSZinv = tr.getVar<int>("cntCSVSZinv");

            float wTT = 1.0;
            float wDY = 1.0;

	    if(cntCSVSZinv == 0)
	    {
		if(njWTTbar_0b)  wTT = 1.0;//njWTTbar_0b->GetBinContent(njWTTbar_0b->FindBin(cntNJetsPt30Eta24Zinv));
		if(njWDYZ_0b)    wDY = njWDYZ_0b->GetBinContent(njWDYZ_0b->FindBin(cntNJetsPt30Eta24Zinv));
	    }
	    else
	    {
		if(njWTTbar_g1b) wTT = 1.0;//njWTTbar_g1b->GetBinContent(njWTTbar_g1b->FindBin(cntNJetsPt30Eta24Zinv));
		if(njWDYZ_g1b)   wDY = njWDYZ_g1b->GetBinContent(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));
	    }
            //std::cout<<wDY<<std::endl;
            float nJet1bfakeWgt = 1.0;
            float nJet2bfakeWgt = 1.0;
            float nJet3bfakeWgt = 1.0;

            if(MCfake1b)   nJet1bfakeWgt = MCfake1b->GetBinContent(MCfake1b->FindBin(cntNJetsPt30Eta24Zinv));
            if(MCfake2b)   nJet2bfakeWgt = MCfake2b->GetBinContent(MCfake2b->FindBin(cntNJetsPt30Eta24Zinv));
            if(MCfake3b)   nJet3bfakeWgt = MCfake3b->GetBinContent(MCfake3b->FindBin(cntNJetsPt30Eta24Zinv));

	    float normWgt0b = ScaleFactors::sf_norm0b();
            float normttbar = ScaleFactorsttBar::sf_norm0b(); 

            tr.registerDerivedVar("nJetWgtTTbar", wTT);
            tr.registerDerivedVar("nJetWgtDYZ",   wDY);

            tr.registerDerivedVar("nJet1bfakeWgt", nJet1bfakeWgt);
            tr.registerDerivedVar("nJet2bfakeWgt", nJet2bfakeWgt);
            tr.registerDerivedVar("nJet3bfakeWgt", nJet3bfakeWgt);

            tr.registerDerivedVar("normWgt0b", normWgt0b);
            tr.registerDerivedVar("normttbar", normttbar);
        }

    public:
        NJetWeight()
        {
            TH1::AddDirectory(false);

            njWTTbar_0b  = nullptr;
            njWDYZ_0b    = nullptr;
            njWTTbar_g1b = nullptr;
            njWDYZ_g1b   = nullptr;
            MCfake1b = nullptr;
            MCfake2b = nullptr;
            MCfake3b = nullptr;

            TFile *f = new TFile("njetWgtHists.root");
            if(f)
            {
                MCfake1b = static_cast<TH1*>(f->Get("h_njRatio_1fake"));
                MCfake2b = static_cast<TH1*>(f->Get("h_njRatio_2fake"));
                MCfake3b = static_cast<TH1*>(f->Get("h_njRatio_3fake"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open: njetWgtHists.root" << std::endl;
            }

            f = new TFile("dataMCweights.root");
            if(f)
            {
                //njWTTbar_0b  = 1;//static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_loose0"));
                njWDYZ_0b    = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_0b_loose0_mt2_MET"));//0b_loose0"));
                //njWTTbar_g1b = 1;//static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_loose0"));
                njWDYZ_g1b   = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_g1b_loose0_mt2_MET"));//g1b_loose0"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open: dataMCweights.root" << std::endl;
            }
        }

        ~NJetWeight()
        {
            //if(njWTTbar) delete njWTTbar;
            //if(njWDYZ)   delete njWDYZ;
        }

        void operator()(NTupleReader& tr)
        {
            generateWeight(tr);
        }
    };

    class LepInfo
    {
    private:
        TRandom3 *tr3;
        void lepInfo(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& muonsLVec    = tr.getVec<TLorentzVector>("muonsLVec");
            const std::vector<float>& muonsRelIso          = tr.getVec<float>("muonsRelIso");
            const std::vector<float>& muonsMiniIso         = tr.getVec<float>("muonsMiniIso");
            const std::vector<float>& muonsCharge          = tr.getVec<float>("muonsCharge");
            const std::vector<TLorentzVector>& jetsLVec     = tr.getVec<TLorentzVector>("jetsLVec");
            const std::vector<float>& recoJetschargedEmEnergyFraction     = tr.getVec<float>("recoJetschargedEmEnergyFraction");
            const std::vector<float>& recoJetschargedHadronEnergyFraction = tr.getVec<float>("recoJetschargedHadronEnergyFraction");
            const std::vector<int> & muonsFlagIDVec = tr.getVec<int>("muonsFlagMedium");
            const std::vector<int>&  elesFlagIDVec  = tr.getVec<int>("elesFlagVeto");

            const std::vector<float>& muonspfActivity      = tr.getVec<float>("muonspfActivity");
            const std::vector<float>& elespfActivity       = tr.getVec<float>("elespfActivity");
            const std::vector<float>& W_emu_pfActivityVec  = tr.getVec<float>("W_emu_pfActivityVec");

            //const float& ht                             = tr.getVar<float>("ht");
            const float& met                            = tr.getVar<float>("met");
            const float& metphi                         = tr.getVar<float>("metphi");

            const std::vector<TLorentzVector, std::allocator<TLorentzVector> > elesLVec = tr.getVec<TLorentzVector>("elesLVec");
            const std::vector<float>& elesMiniIso          = tr.getVec<float>("elesMiniIso");
            const std::vector<float>& elesCharge           = tr.getVec<float>("elesCharge");
            const std::vector<unsigned int>& elesisEB       = tr.getVec<unsigned int>("elesisEB",true);

            //const int& nTopCandSortedCnt = tr.getVar<int>("nTopCandSortedCntZinv");

            bool passMuonVeto = false;
            bool passEleVeto = false;

            if(tr.checkBranch("passMuonVeto"))
            {
                passMuonVeto = tr.getVar<bool>("passMuonVeto");
            }
            if(tr.checkBranch("passEleVeto"))
            {
                passEleVeto = tr.getVar<bool>("passEleVeto");
            }

            std::vector<const TLorentzVector*>* genMatchIsoElecInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genMatchElecInAcc = new std::vector<const TLorentzVector*>();
            std::vector<float>* genMatchElecInAccRes = new std::vector<float>();
            std::vector<const TLorentzVector*>* genElecInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genElec = new std::vector<const TLorentzVector*>();
            std::vector<float>* genMatchIsoElecInAccAct = new std::vector<float>();
            std::vector<float>* genMatchElecInAccAct = new std::vector<float>();
            std::vector<float>* genElecInAccAct = new std::vector<float>();
            std::vector<float>* genElecAct = new std::vector<float>();

            std::vector<TLorentzVector>* cutElecVec = new std::vector<TLorentzVector>();
            std::vector<float>* cutElecCharge = new std::vector<float>();
            std::vector<float>* cutElecActivity = new std::vector<float>();

            std::vector<TLorentzVector> cutElecVecRecoOnly;

            std::vector<const TLorentzVector*>* genMatchIsoMuInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genMatchMuInAcc = new std::vector<const TLorentzVector*>();
            std::vector<float>* genMatchMuInAccRes = new std::vector<float>();
            std::vector<const TLorentzVector*>* genMuInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genMu = new std::vector<const TLorentzVector*>();
            std::vector<float>* genMatchIsoMuInAccAct = new std::vector<float>();
            std::vector<float>* genMatchMuInAccAct = new std::vector<float>();
            std::vector<float>* genMuInAccAct = new std::vector<float>();
            std::vector<float>* genMuAct = new std::vector<float>();
            std::vector<TLorentzVector>* cutMuVec = new std::vector<TLorentzVector>();
            std::vector<float>* cutMuCharge = new std::vector<float>();
            std::vector<float>* cutMuActivity = new std::vector<float>();

            //std::vector<TLorentzVector>* Zrecopt = new std::vector<TLorentzVector>();

            std::vector<TLorentzVector> cutMuVecRecoOnly;

            //muon selections
            int sumMuCharge = 0;
            int nTriggerMuons = 0;

            for(int i = 0; i < muonsLVec.size(); ++i)
            {
                if(AnaFunctions::passMuon( muonsLVec[i], 0.0, 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr)) // emulates muons with pt but no iso requirements.
                {
                    cutMuVecRecoOnly.push_back(muonsLVec[i]);
                }
                if(AnaFunctions::passMuon( muonsLVec[i], muonsMiniIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                {
                    if(AnaFunctions::passMuon( muonsLVec[i], muonsRelIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                    {
                        if(nTriggerMuons == 0 && muonsLVec[i].Pt() > 17)  nTriggerMuons++;
                        else if(muonsLVec[i].Pt() > 8)  nTriggerMuons++;
                    }
                    cutMuVec->push_back(muonsLVec[i]);
                    //std::cout<<"cutMuVec PT "<<muonsLVec[i].Pt()<<std::endl; 
                    cutMuCharge->push_back(muonsCharge[i]);
                    cutMuActivity->push_back(muonspfActivity[i]);
                    if(muonsCharge[i] > 0) sumMuCharge++;
                    else                   sumMuCharge--;
                }
            }
            //std::cout<<"New Muon Selection "<<(*cutMuVec).size()<<std::endl;
            //electron selection
            int sumElecCharge = 0;
            for(int i = 0; i < elesLVec.size(); ++i)
            {
                if(AnaFunctions::passElectron(elesLVec[i], 0.0, -1, elesisEB[i], elesFlagIDVec[i], AnaConsts::elesMiniIsoArr)) // emulates muons with pt but no iso requirements.
                {
                    cutElecVecRecoOnly.push_back(elesLVec[i]);
                }

                if(AnaFunctions::passElectron(elesLVec[i], elesMiniIso[i], -1, elesisEB[i], elesFlagIDVec[i], AnaConsts::elesMiniIsoArr))
                {
                    cutElecVec->push_back(elesLVec[i]);
                    cutElecCharge->push_back(elesCharge[i]);
                    cutElecActivity->push_back(AnaFunctions::getElectronActivity(elesLVec[i], jetsLVec, recoJetschargedHadronEnergyFraction, AnaConsts::elesAct));
                    if(elesCharge[i] > 0) sumElecCharge++;
                    else                  sumElecCharge--;
                }
            }


            //mu45 non-iso trigger emulation
            const float effsnom2012ABC[] = {0.928,0.8302,0.8018};
            const float upedge2012ABC[] = { 0.9, 1.2, 2.1};
            bool muTrigMu45 = false;
            for(TLorentzVector& mu : *cutMuVec)
            {
                if(mu.Pt() > 50)
                {
                    for(int iBin = 0; iBin < sizeof(effsnom2012ABC)/sizeof(float); ++iBin)
                    {
                        if(mu.Eta() < upedge2012ABC[iBin] && tr3 && tr3->Uniform() < effsnom2012ABC[iBin])
                        {
                            muTrigMu45 = true;
                            break;
                        }
                    }
                    if(muTrigMu45) break;
                }
            }

            float genHt = 0.0;

            const float   minMuPt = 20.0,   highMuPt = 50.0;
            const float minElecPt = 33.0, highElecPt = 33.0;
            float nuPt1 = -999.9, nuPt2 = -999.9;

            //Gen info parsing
            float genZPt = -999.9, genZEta = -999.9, genZmass = -999.9, genZPhi;
            int nZ = 0;
            TLorentzVector genZ;
            int pdgIdZDec = 0;
            if(tr.checkBranch("genDecayLVec") && tr.checkBranch("genDecayPdgIdVec") && tr.checkBranch("genDecayPdgIdVec") && tr.checkBranch("W_emuVec"))
            {
                const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
                const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");
                const std::vector<int>& W_emuVec                = tr.getVec<int>("W_emuVec");

                for(int i = 0; i < genDecayPdgIdVec.size() && i < genDecayLVec.size(); ++i)
                {
                    if((abs(genDecayPdgIdVec[i]) != 0 &&  abs(genDecayPdgIdVec[i]) < 6) || (abs(genDecayPdgIdVec[i]) > 100 && abs(genDecayPdgIdVec[i]) < 10000)) genHt += genDecayLVec[i].Pt();

                    if(genDecayPdgIdVec[i] ==  13) nuPt1 = genDecayLVec[i].Pt();
                    if(genDecayPdgIdVec[i] == -13) nuPt2 = genDecayLVec[i].Pt();
                }

                for(int index = 0; index < W_emuVec.size(); ++index)
                {
                    int i = W_emuVec[index];

                    //muon efficiency and acceptance
                    if(abs(genDecayPdgIdVec[i]) == 13)
                    {
                        genMu->push_back(&genDecayLVec[i]);
                        genMuAct->push_back(W_emu_pfActivityVec[index]);
                        if(AnaFunctions::passMuonAccOnly(genDecayLVec[i], AnaConsts::muonsMiniIsoArr) && genDecayLVec[i].Pt() > minMuPt)
                        {
                            genMuInAcc->push_back(&genDecayLVec[i]);
                            genMuInAccAct->push_back(genMuAct->back());
                            float dRMin = 999.9;
                            float matchPt = -999.9;
                            for(int j = 0; j < cutMuVecRecoOnly.size(); ++j)
                            {
                                float dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutMuVecRecoOnly[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                    matchPt = cutMuVecRecoOnly[j].Pt();
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchMuInAcc->push_back(&genDecayLVec[i]);
                                genMatchMuInAccAct->push_back(genMuAct->back());
                                genMatchMuInAccRes->push_back((genDecayLVec[i].Pt() - matchPt)/genDecayLVec[i].Pt());
                            }

                            dRMin = 999.9;
                            for(int j = 0; j < cutMuVec->size(); ++j)
                            {
                                float dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], (*cutMuVec)[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchIsoMuInAcc->push_back(&genDecayLVec[i]);
                                genMatchIsoMuInAccAct->push_back(genMuAct->back());
                            }
                        }
                    }

                    //Elec efficiency and acceptance
                    if(abs(genDecayPdgIdVec[i]) == 11)
                    {
                        genElec->push_back(&genDecayLVec[i]);
                        genElecAct->push_back(W_emu_pfActivityVec[index]);
                        if(AnaFunctions::passElectronAccOnly(genDecayLVec[i], AnaConsts::elesMiniIsoArr) && genDecayLVec[i].Pt() > 20)
                        {
                            genElecInAcc->push_back(&genDecayLVec[i]);
                            genElecInAccAct->push_back(genElecAct->back());
                            float dRMin = 999.9;
                            float matchPt = -999.9;
                            for(int j = 0; j < cutElecVecRecoOnly.size(); ++j)
                            {
                                float dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutElecVecRecoOnly[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                    matchPt = cutElecVecRecoOnly[j].Pt();
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchElecInAcc->push_back(&genDecayLVec[i]);
                                genMatchElecInAccAct->push_back(genElecAct->back());
                                genMatchElecInAccRes->push_back((genDecayLVec[i].Pt() - matchPt)/genDecayLVec[i].Pt());
                            }

                            dRMin = 999.9;
                            for(int j = 0; j < cutElecVec->size(); ++j)
                            {
                                float dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], (*cutElecVec)[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchIsoElecInAcc->push_back(&genDecayLVec[i]);
                                genMatchIsoElecInAccAct->push_back(genElecAct->back());
                            }
                        }
                    }
                }

                for(int j = 0; j <  genDecayPdgIdVec.size(); ++j)
                {
                    if(abs(genDecayPdgIdVec[j]) == 23)
                    {
                        nZ++;
                        genZ = genDecayLVec[j];
                        genZPt = genDecayLVec[j].Pt();
                        genZEta = genDecayLVec[j].Eta();
                        genZPhi = genDecayLVec[j].Phi();
                        genZmass = genDecayLVec[j].M();
                    }
                }
                if(nZ > 1) std::cout << "!!!WARNING MORE THAN 1 Z FOUND!!!" << std::endl;

                if(W_emuVec.size() == 0) pdgIdZDec = 15;
                else if(W_emuVec.size() == 2)
                {
                    if(abs(genDecayPdgIdVec[W_emuVec[0]]) == 11) pdgIdZDec = 11;
                    else if(abs(genDecayPdgIdVec[W_emuVec[0]]) == 13) pdgIdZDec = 13;
                }

            }

            bool passDiMuTrig  = nTriggerMuons >= 2;

            const float zMassMin = 81.0;
            const float zMass    = 91.0;
            const float zMassMax = 101.0;

            float zMuMassCurrent = 1.0e300, zEff = 1.0e100, zAcc = 1.0e100;
            TLorentzVector bestRecoMuZ;
            for(int i = 0; i < cutMuVec->size(); ++i)
            {
                if((*cutMuVec)[i].Pt() < minMuPt) continue;
                for(int j = 0; j < i && j < cutMuVec->size(); ++j)
                {
                    if((*cutMuVec)[j].Pt() < minMuPt) continue;
                    float zm = ((*cutMuVec)[i] + (*cutMuVec)[j]).M();
                    if(fabs(zm - zMass) < fabs(zMuMassCurrent - zMass))
                    {
                        bestRecoMuZ = (*cutMuVec)[i] + (*cutMuVec)[j];
                        zMuMassCurrent = zm;
                    }
                }
            }

            float zElecMassCurrent = 1.0e300;
            TLorentzVector bestRecoElecZ;
            for(int i = 0; i < cutElecVec->size(); ++i)
            {
                if((*cutElecVec)[i].Pt() < minElecPt) continue;
                for(int j = 0; j < i && j < cutElecVec->size(); ++j)
                {
                    if((*cutElecVec)[j].Pt() < minElecPt) continue;
                    float zm = ((*cutElecVec)[i] + (*cutElecVec)[j]).M();
                    //if(zm > zMassMin && zm < zMassMax && fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    if(fabs(zm - zMass) < fabs(zElecMassCurrent - zMass))
                    {
                        bestRecoElecZ = (*cutElecVec)[i] + (*cutElecVec)[j];
                        zElecMassCurrent = zm;
                    }
                }
            }

            float zElMuMassCurrent = 1.0e300;
            TLorentzVector bestRecoElMuZ;
            for(int i = 0; i < cutMuVec->size(); ++i)
            {
                if((*cutMuVec)[i].Pt() < minMuPt) continue;
                for(int j = 0; j < cutElecVec->size(); ++j)
                {
                    if((*cutElecVec)[j].Pt() < minMuPt) continue;
                    float zm = ((*cutMuVec)[i] + (*cutElecVec)[j]).M();
                    //if(zm > zMassMin && zm < zMassMax && fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    if(fabs(zm - zMass) < fabs(zElMuMassCurrent - zMass))
                    {
                        bestRecoElMuZ = (*cutMuVec)[i] + (*cutElecVec)[j];
                        zElMuMassCurrent = zm;
                    }
                }
            }

            TLorentzVector metV, metZ;
            metV.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

            TLorentzVector bestRecoZ = (true/*fabs(bestRecoElecZ.M() - zMass) > fabs(bestRecoMuZ.M() - zMass)*/)?(bestRecoMuZ):(bestRecoElecZ);
            //if(fabs(bestRecoZ.M() - zMass) > fabs(bestRecoElMuZ.M() - zMass)) bestRecoZ = bestRecoElMuZ;

            metZ.SetPtEtaPhiM(bestRecoZ.Pt(), 0.0, bestRecoZ.Phi(), 0.0);
            TLorentzVector cleanMet = metV + metZ;
            //std::cout<<"metZ "<<metZ.Pt()<<std::endl;
            //std::cout<<"metV "<<metV.Pt()<<std::endl;
            bool passDiMuSel   =  passEleVeto && (cutMuVec->size() == 2   && sumMuCharge == 0   && (*cutMuVec)[0].Pt() > highMuPt     && (*cutMuVec)[1].Pt() > minMuPt);
            bool passDiElecSel = passMuonVeto && (cutElecVec->size() == 2 && sumElecCharge == 0 && (*cutElecVec)[0].Pt() > highElecPt && (*cutElecVec)[1].Pt() > minElecPt);
            bool passElMuSel = (cutMuVec->size() == 1 && cutElecVec->size() == 1 && sumElecCharge == -sumMuCharge && (*cutMuVec)[0].Pt() > highMuPt && (*cutElecVec)[0].Pt() > minMuPt);

            bool passMuZinvSel   =  passEleVeto && (cutMuVec->size() == 2   && sumMuCharge == 0   && (*cutMuVec)[0].Pt() > highMuPt     && (*cutMuVec)[1].Pt() > minMuPt)     && (bestRecoMuZ.M() > zMassMin)   && (bestRecoMuZ.M() < zMassMax);
            bool passElecZinvSel = passMuonVeto && (cutElecVec->size() == 2 && sumElecCharge == 0 && (*cutElecVec)[0].Pt() > highElecPt && (*cutElecVec)[1].Pt() > minElecPt) && (bestRecoElecZ.M() > zMassMin) && (bestRecoElecZ.M() < zMassMax);
            bool passElMuZinvSel = (cutMuVec->size() == 1 && cutElecVec->size() == 1 && sumElecCharge == -sumMuCharge && (*cutMuVec)[0].Pt() > highMuPt && (*cutElecVec)[0].Pt() > minMuPt) && (bestRecoElMuZ.M() > zMassMin) && (bestRecoElMuZ.M() < zMassMax);

            float cutMuPt1 = -999.9;
            float cutMuPt2 = -999.9;
            if(cutMuVec->size() >= 1) cutMuPt1 = cutMuVec->at(0).Pt();
            if(cutMuVec->size() >= 2) cutMuPt2 = cutMuVec->at(1).Pt();

            float cutElecPt1 = -999.9;
            float cutElecPt2 = -999.9;
            if(cutElecVec->size() >= 1) cutElecPt1 = cutElecVec->at(0).Pt();
            if(cutElecVec->size() >= 2) cutElecPt2 = cutElecVec->at(1).Pt();

            float mindPhiMetJ = 999.9;
            int jc = 0;
            for(const TLorentzVector& jet : jetsLVec)
            {
                if(jc >= 3) break;
                jc++;
                mindPhiMetJ = std::min(mindPhiMetJ, (float)fabs(ROOT::Math::VectorUtil::DeltaPhi(genZ, jet)));
            }

            float bestRecoZPt = bestRecoZ.Pt();
            float cleanMetPt = cleanMet.Pt();

            tr.registerDerivedVar("bestRecoZPt", bestRecoZPt);
            tr.registerDerivedVar("bestRecoZM", bestRecoZ.M());
            tr.registerDerivedVar("cleanMetPt", cleanMetPt);
            tr.registerDerivedVar("cleanMetPhi", cleanMet.Phi());
            //tr.registerDerivedVar("cleanMet2Pt", cleanMet2Pt);
            tr.registerDerivedVar("genHt", genHt);
            tr.registerDerivedVar("cutMuPt1", cutMuPt1);
            tr.registerDerivedVar("cutMuPt2", cutMuPt2);
            tr.registerDerivedVar("cutElecPt1", cutElecPt1);
            tr.registerDerivedVar("cutElecPt2", cutElecPt2);
            tr.registerDerivedVar("mindPhiMetJ", mindPhiMetJ);

            tr.registerDerivedVar("ZPtRes", (bestRecoZPt - genZPt)/genZPt);
            tr.registerDerivedVar("ZEtaRes", bestRecoZ.Eta() - genZEta);
            tr.registerDerivedVar("ZPhiRes", bestRecoZ.Phi() - genZPhi);
            tr.registerDerivedVar("ZMRes", (bestRecoZ.M() - genZmass)/genZmass);

            tr.registerDerivedVec("cutMuVec", cutMuVec);
            tr.registerDerivedVec("cutElecVec", cutElecVec);
            tr.registerDerivedVec("cutMuActivity", cutMuActivity);
            tr.registerDerivedVec("cutElecActivity", cutElecActivity);
            tr.registerDerivedVec("genMu", genMu);
            tr.registerDerivedVar("ngenMu", static_cast<float>(genMu->size()));
            tr.registerDerivedVec("genMuInAcc", genMuInAcc);
            tr.registerDerivedVec("genMuAct", genMuAct);
            tr.registerDerivedVar("ngenMuInAcc", static_cast<float>(genMuInAcc->size()));
            tr.registerDerivedVec("genMuInAccAct", genMuInAccAct);
            tr.registerDerivedVec("genMatchMuInAcc", genMatchMuInAcc);
            tr.registerDerivedVec("genMatchMuInAccRes", genMatchMuInAccRes);
            tr.registerDerivedVec("genMatchIsoMuInAcc", genMatchIsoMuInAcc);
            tr.registerDerivedVar("ngenMatchMuInAcc", static_cast<float>(genMatchMuInAcc->size()));
            tr.registerDerivedVec("genMatchMuInAccAct", genMatchMuInAccAct);
            tr.registerDerivedVec("genMatchIsoMuInAccAct", genMatchIsoMuInAccAct);

            tr.registerDerivedVec("genElec", genElec);
            tr.registerDerivedVar("ngenElec", static_cast<float>(genElec->size()));
            tr.registerDerivedVec("genElecInAcc", genElecInAcc);
            tr.registerDerivedVec("genElecAct", genElecAct);
            tr.registerDerivedVar("ngenElecInAcc", static_cast<float>(genElecInAcc->size()));
            tr.registerDerivedVec("genElecInAccAct", genElecInAccAct);
            tr.registerDerivedVec("genMatchElecInAcc", genMatchElecInAcc);
            tr.registerDerivedVec("genMatchElecInAccRes", genMatchElecInAccRes);
            tr.registerDerivedVec("genMatchIsoElecInAcc", genMatchIsoElecInAcc);
            tr.registerDerivedVar("ngenMatchElecInAcc", static_cast<float>(genMatchElecInAcc->size()));
            tr.registerDerivedVec("genMatchElecInAccAct", genMatchElecInAccAct);
            tr.registerDerivedVec("genMatchIsoElecInAccAct", genMatchIsoElecInAccAct);

            tr.registerDerivedVar("genZPt", genZPt);
            tr.registerDerivedVar("genZEta", genZEta);
            tr.registerDerivedVar("genZmass", genZmass);
            tr.registerDerivedVar("pdgIdZDec", pdgIdZDec);
            tr.registerDerivedVar("passDiMuIsoTrig", passDiMuTrig);
            tr.registerDerivedVar("passSingleMu45", muTrigMu45);

            tr.registerDerivedVar("passDiMuSel", passDiMuSel);
            tr.registerDerivedVar("passDiElecSel", passDiElecSel);
	    tr.registerDerivedVar("passElMuSel", passElMuSel);

            tr.registerDerivedVar("passMuZinvSel", passMuZinvSel);
            tr.registerDerivedVar("passElecZinvSel", passElecZinvSel);
            tr.registerDerivedVar("passElMuZinvSel", passElMuZinvSel);

            tr.registerDerivedVar("Zrecopt",bestRecoZPt);
        }

    public:
        LepInfo()
        {
            tr3 = new TRandom3();
        }

        void operator()(NTupleReader& tr)
        {
            lepInfo(tr);
        }
    };

    class Fakebtagvectors
    {
    private:
        void fakebtagvectors(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& jetsLVecLepCleaned = tr.getVec<TLorentzVector>("jetsLVecLepCleaned");
            //New code
            const std::vector<float>& cleanJetpt30ArrBTag = tr.getVec<float>("recoJetsBtag_0_LepCleaned");
//            const std::vector<float>& cleanJetpt30ArrBTag = tr.getVec<float>("recoJetsBtag_0_LepCleaned");

            float maxCSV = 0.0;
            float secCSV = 0.0;
            float tenCSV = 0.0;
            int iMaxCSV = -1;
            int iSecCSV = -1;
            int iTenCSV = -1;

            if(jetsLVecLepCleaned.size() != cleanJetpt30ArrBTag.size()) std::cout << "fakebtagvectors(...): Vector size missmatch!!!!" << std::endl;

            int njet = 0;

            //find index of 3 highest CSV values
            for(int i = 0; i < cleanJetpt30ArrBTag.size(); ++i)
            {
                //Skip jets which cannot pass bTag Acceptance requirements
                if(!AnaFunctions::jetPassCuts(jetsLVecLepCleaned[i], AnaConsts::bTagArr)) continue;

                //count possible fake b-jets
                njet++;

                if(cleanJetpt30ArrBTag[i] > maxCSV)
                {
                    tenCSV = secCSV;
                    secCSV = maxCSV;
                    maxCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = iSecCSV;
                    iSecCSV = iMaxCSV;
                    iMaxCSV = i;
                }
                else if(cleanJetpt30ArrBTag[i] > secCSV)
                {
                    tenCSV = secCSV;
                    secCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = iSecCSV;
                    iSecCSV = i;
                }
                else if(cleanJetpt30ArrBTag[i] > tenCSV)
                {
                    tenCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = i;
                }
            }

            std::vector<float>* cleanJetpt30ArrBTag1fake = new std::vector<float>(cleanJetpt30ArrBTag);
            std::vector<float>* cleanJetpt30ArrBTag2fake = new std::vector<float>(cleanJetpt30ArrBTag);
            std::vector<float>* cleanJetpt30ArrBTag3fake = new std::vector<float>(cleanJetpt30ArrBTag);
            std::vector<float>* fakedCSVValues = new std::vector<float>();

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag1fake)[iMaxCSV] = 0.99;

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag2fake)[iMaxCSV] = 0.99;
            if(iSecCSV >= 0) (*cleanJetpt30ArrBTag2fake)[iSecCSV] = 0.99;

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iMaxCSV] = 0.99;
            if(iSecCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iSecCSV] = 0.99;
            if(iTenCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iTenCSV] = 0.99;

            if(iMaxCSV >= 0) fakedCSVValues->push_back(maxCSV);
            if(iSecCSV >= 0) fakedCSVValues->push_back(secCSV);
            if(iTenCSV >= 0) fakedCSVValues->push_back(tenCSV);

            //Calculate the combinatoric weights for b-jet faking
            float weight1fakeb = TMath::Binomial(njet, 1);
            float weight2fakeb = TMath::Binomial(njet, 2);
            float weight3fakeb = TMath::Binomial(njet, 3);
            //check for nans
            if(weight1fakeb != weight1fakeb) weight1fakeb = 0.0;
            if(weight2fakeb != weight2fakeb) weight2fakeb = 0.0;
            if(weight3fakeb != weight3fakeb) weight3fakeb = 0.0;

            tr.registerDerivedVar("weight1fakeb", weight1fakeb);
            tr.registerDerivedVar("weight2fakeb", weight2fakeb);
            tr.registerDerivedVar("weight3fakeb", weight3fakeb);

            tr.registerDerivedVec("cleanJetpt30ArrBTag1fake", cleanJetpt30ArrBTag1fake);
            tr.registerDerivedVec("cleanJetpt30ArrBTag2fake", cleanJetpt30ArrBTag2fake);
            tr.registerDerivedVec("cleanJetpt30ArrBTag3fake", cleanJetpt30ArrBTag3fake);
            tr.registerDerivedVec("fakedCSVValues", fakedCSVValues);
            tr.registerDerivedVar("maxCSV", maxCSV);
        }
    public:

        void operator()(NTupleReader& tr)
        {
            fakebtagvectors(tr);
        }
    };

    class GetSearchBin
    {
    private:
        SearchBins sbins;

        void getSearchBin(NTupleReader& tr)
        {
            const int& cntCSVS = tr.getVar<int>("cntCSVSZinv");
            const int& nTopCandSortedCnt = tr.getVar<int>("nTopCandSortedCntZinv");
            //const int& nTopCandSortedCnt1b = tr.getVar<int>("nTopCandSortedCntZinv1b");
            //const int& nTopCandSortedCnt2b = tr.getVar<int>("nTopCandSortedCntZinv2b");
            //const int& nTopCandSortedCnt3b = tr.getVar<int>("nTopCandSortedCntZinv3b");
            const float& cleanMet = tr.getVar<float>("cleanMetPt");
            const float& cleanMetPhi = tr.getVar<float>("cleanMetPhi");
            const float& MT2 = tr.getVar<float>("best_had_brJet_MT2Zinv");
            //const float& MT2_1b = tr.getVar<float>("best_had_brJet_MT2Zinv1b");
            //const float& MT2_2b = tr.getVar<float>("best_had_brJet_MT2Zinv2b");
            //const float& MT2_3b = tr.getVar<float>("best_had_brJet_MT2Zinv3b");
            //const float& weight1fakeb = tr.getVar<float>("weight1fakeb");
            //const float& weight2fakeb = tr.getVar<float>("weight2fakeb");
            //const float& weight3fakeb = tr.getVar<float>("weight3fakeb");
            //
            //const float& nJet1bfakeWgt = tr.getVar<float>("nJet1bfakeWgt");
            //const float& nJet2bfakeWgt = tr.getVar<float>("nJet2bfakeWgt");
            //const float& nJet3bfakeWgt = tr.getVar<float>("nJet3bfakeWgt");
            const float& HT            = tr.getVar<float>("HTZinv");
            //const int& nTopCandSortedCnt = tr.getVar<int>("nTopCandSortedCntZinv");
            //            //top
            /*
            std::shared_ptr<TopTagger> ttPtr;
            int monoJet;
           const TopTaggerResults& ttr = ttPtr->getResults();
               std::vector<TopObject*> Ntop = ttr.getTops();
              for(int i=1; i<nTopCandSortedCnt; i++){
               if(Ntop[i]->getNConstituents() == 1) monoJet++;
                  }
             std::cout<<monoJet<<std::endl;
            */
            //int nSearchBin = sbins.find_Binning_Index(cntCSVS, nTopCandSortedCnt, MT2, cleanMet);
            int nSearchBin = sbins.find_Binning_Index(cntCSVS, nTopCandSortedCnt, MT2, cleanMet, HT);            

            //std::vector<std::pair<float, float> > * nb0Bins = new std::vector<std::pair<float, float> >();
            //std::vector<std::pair<float, float> > * nb0NJwBins = new std::vector<std::pair<float, float> >();
            //std::vector<float> * nb0BinsNW = new std::vector<float>();

            //weights based on total N(b) yields vs. N(b) = 0 control region
            //These weights are derived from the rato of events in the N(t) = 1, 2, 3 bins after all baseline cuts except b tag between the
            //N(b) = 0 control region and each N(b) signal region using Z->nunu MC.  They account for both the combinatoric reweighting factor
            //as well as the different event yields between the control region and each signal region.
            //const float wnb01 = 3.820752e-02;//3.478840e-02;//6.26687e-2;
            //const float wnb02 = 2.946461e-03;//2.586369e-03;//5.78052e-3;
            //const float wnb03 = 1.474770e-04;//1.640077e-04;//7.08235e-4;
            //
            //// weights to apply when doing b-faking
            //const float w1b = wnb01 * weight1fakeb;
            //const float w2b = wnb02 * weight2fakeb;
            //const float w3b = wnb03 * weight3fakeb;

            if(cntCSVS == 0)
            {
                //nb0Bins->push_back(std::make_pair(find_Binning_Index(0, nTopCandSortedCnt, MT2, cleanMet), 1.0));
                //nb0Bins->emplace_back(std::pair<float, float>(sbins.find_Binning_Index(1, nTopCandSortedCnt1b, MT2_1b, cleanMet, HT), wnb01 * weight1fakeb));
                //nb0Bins->emplace_back(std::pair<float, float>(sbins.find_Binning_Index(2, nTopCandSortedCnt2b, MT2_2b, cleanMet, HT), wnb02 * weight2fakeb));
                //nb0Bins->emplace_back(std::pair<float, float>(sbins.find_Binning_Index(3, nTopCandSortedCnt3b, MT2_3b, cleanMet, HT), wnb03 * weight3fakeb));

                //nb0NJwBins->emplace_back(std::pair<float, float>(sbins.find_Binning_Index(1, nTopCandSortedCnt1b, MT2_1b, cleanMet, HT), nJet1bfakeWgt));
                //nb0NJwBins->emplace_back(std::pair<float, float>(sbins.find_Binning_Index(2, nTopCandSortedCnt2b, MT2_2b, cleanMet, HT), nJet2bfakeWgt));
                //nb0NJwBins->emplace_back(std::pair<float, float>(sbins.find_Binning_Index(3, nTopCandSortedCnt3b, MT2_3b, cleanMet, HT), nJet3bfakeWgt));

                //nb0BinsNW->emplace_back(sbins.find_Binning_Index(1, nTopCandSortedCnt1b, MT2_1b, cleanMet, HT));
                //nb0BinsNW->emplace_back(sbins.find_Binning_Index(2, nTopCandSortedCnt2b, MT2_2b, cleanMet, HT));
                //nb0BinsNW->emplace_back(sbins.find_Binning_Index(3, nTopCandSortedCnt3b, MT2_3b, cleanMet, HT));
            }

            tr.registerDerivedVar("nSearchBin", nSearchBin);
            //tr.registerDerivedVec("nb0BinsNW", nb0BinsNW);
            //tr.registerDerivedVec("nb0Bins", nb0Bins);
            //tr.registerDerivedVec("nb0NJwBins", nb0NJwBins);
            //tr.registerDerivedVar("weight1fakebComb", w1b);
            //tr.registerDerivedVar("weight2fakebComb", w2b);
            //tr.registerDerivedVar("weight3fakebComb", w3b);
        }

    public:

        GetSearchBin(std::string sb_era) : sbins(sb_era) {}

        void operator()(NTupleReader& tr)
        {
            getSearchBin(tr);
        }
    };


    class TriggerInfo
    {
    private:
	int indexMuTrigger;
	int indexElecTrigger;
        int indexMETMHTTrigger;
        bool miniTuple_, noMC_;

	float GetMuonTriggerEff(const float& muEta) 
	{
            if     (-2.6 <= muEta && muEta < -2.2) return 0.8020833;
            else if(-2.2 <= muEta && muEta < -1.8) return 0.8113949;
            else if(-1.8 <= muEta && muEta < -1.4) return 0.8111837;
            else if(-1.4 <= muEta && muEta < -1.0) return 0.8824405;
            else if(-1.0 <= muEta && muEta < -0.6) return 0.9024091;
            else if(-0.6 <= muEta && muEta < -0.2) return 0.8737864;
            else if(-0.2 <= muEta && muEta <  0.2) return 0.9186085;
            else if( 0.2 <= muEta && muEta <  0.6) return 0.8759649;
            else if( 0.6 <= muEta && muEta <  1.0) return 0.8940410;
            else if( 1.0 <= muEta && muEta <  1.4) return 0.8848286;
            else if( 1.4 <= muEta && muEta <  1.8) return 0.8293217;
            else if( 1.8 <= muEta && muEta <  2.2) return 0.8263979;
            else if( 2.2 <= muEta && muEta <  2.6) return 0.7605634;
            else                                   return 0.000;
	}

	float GetTriggerEffWeight(const float& met, const float& ht) 
	{
	    if (ht<1000)
	    {
		if (met<25) return 0.001542561;
		else if (met<50) return 0.003222389;
		else if (met<75) return 0.00987073;
		else if (met<100) return 0.03865682;
		else if (met<125) return 0.1387231;
		else if (met<150) return 0.3564816;
		else if (met<175) return 0.6276442;
		else if (met<200) return 0.8154821;
		else if (met<275) return 0.9340538;
		else if (met<400) return 0.9858562; 
                else if (met<600) return 0.9931507;
                else if (met<1000) return 1.00;
		else return 1.00;
	    } 
	    else 
	    {
		if (met<25) return  0.02067183;
		else if (met<50) return 0.02504944;
		else if (met<75) return 0.04486466;
		else if (met<100) return 0.07434402;
		else if (met<125) return 0.1518288;
		else if (met<150) return 0.2802669;
		else if (met<175) return 0.4642409;
		else if (met<200) return 0.6596434;
		else if (met<275) return 0.8510453;
		else if (met<400) return 0.9563492;
                else if (met<600) return 0.9874214;
                else if (met<1000) return 0.9736842; 
		else return 0.9736842;
	    }
	}
	float GetTriggerEffStatUncHi(const float& met, const float& ht) 
	{
	    if (ht<1000)
	    {
		if (met<25) return 0.0001251554;
		else if (met<50) return 0.0001310897;
		else if (met<75) return 0.0002597269;
		else if (met<100) return 0.0006525702;
		else if (met<125) return 0.001545856;
		else if (met<150) return 0.002821274;
                else if (met<200) return 0.003691577;
		else if (met<275) return 0.003877182;
		else if (met<400) return 0.002294442; 
                else if (met<600) return 0.002045071;
                else if (met<1000) return 0.003725375;
		else return 0.00;
	    } 
	    else 
	    {
		if (met<25) return 0.004283915;
		else if (met<50) return 0.003169914;
		else if (met<75) return 0.004349597;
		else if (met<100) return 0.006241982;
		else if (met<125) return 0.01001983;
		else if (met<150) return 0.01455422;
		else if (met<175) return 0.0183275;
		else if (met<200) return 0.01960093;
		else if (met<275) return 0.01062354;
		else if (met<400) return 0.007445741;
                else if (met<600) return 0.006010458;
                else if (met<1000) return 0.01697945; 
		else return 0.01697945;
	    }
	}
	float GetTriggerEffStatUncLo(const float& met, const float& ht) 
	{
	    if (ht<1000)
	    {
                if (met<25) return 0.0001160878;
                else if (met<50) return 0.000126075;
                else if (met<75) return 0.0002532144;
                else if (met<100) return 0.000642253;
                else if (met<125) return 0.001531628;
                else if (met<150) return 0.002811409;
                else if (met<175) return 0.003706407;
                else if (met<200) return 0.003940439;
                else if (met<275) return 0.00236968;
                else if (met<400) return 0.002358961;
                else if (met<600) return 0.006617554;
                else if (met<1000) return 0.1422293;
                else return 0.1422293;
            }
	   
	    else 
	    {
		if (met<25) return 0.003609465;
		else if (met<50) return 0.002838673;
		else if (met<75) return 0.003996443;
		else if (met<100) return 0.005811049;
		else if (met<125) return 0.009521872;
		else if (met<150) return 0.01412113;
		else if (met<175) return 0.01823465;
		else if (met<200) return 0.02013986;
		else if (met<275) return 0.01126014;
		else if (met<400) return 0.008759573;
                else if (met<600) return 0.009833846;
                else if (met<1000) return 0.03365661; 
		else return 0.03365661;
	    }
	}
	float GetTriggerEffSystUncHi(const float& met, const float& ht) 
	{
	    return 0.0;
	    /* if (met<100) return 0.0272; */
	    /* else if (met<150) return 0.0872; */
	    /* else if (met<175) return 0.1505; */
	    /* else if (met<200) return 0.0423; */
	    /* else if (met<275) return 0.0112; */
	    /* else if (met<400) return 0.0001;  */
	    /* else return 0.0; */
	}
	float GetTriggerEffSystUncLo(const float& met, const float& ht) 
	{
	    return 0.0;
	    /* if (met<100) return 0.0120; */
	    /* else if (met<150) return 0.0872; */
	    /* else if (met<175) return 0.1505; */
	    /* else if (met<200) return 0.0792; */
	    /* else if (met<275) return 0.0112; */
	    /* else if (met<400) return 0.0001;  */
	    /* else return 0.0018; */
	}

        void triggerInfo(NTupleReader& tr)
        {
            const std::vector<std::string>& triggerNames = tr.getVec<std::string>("TriggerNames");
            const std::vector<int>& passTrigger          = tr.getVec<int>("PassTrigger");

            bool passMuTrigger = false;
            bool passElecTrigger = false;
            bool passMETMHTTrigger = false;

	    const std::string muTrigName = "HLT_Mu50_v";//"HLT_Mu45_eta2p1_v";
	    const std::string elecTrigName = "HLT_floatEle33_CaloIdL_GsfTrkIdVL_MW_v";
            const std::string metmhtTrigName = "HLT_PFMET110_PFMHT110_IDTight_v";

            // Find the index of our triggers if we don't know them already
            //if(indexMuTrigger == -1 || indexElecTrigger == -1 || indexMETMHTTrigger == -1)
            //{
            //    for(int i = 0; i < triggerNames.size(); ++i)
            //    {
            //        //if(triggerNames[i].find(muTrigName) != std::string::npos)
            //        //{
            //        //    indexMuTrigger = i;
            //        //}
            //        if(triggerNames[i].find(elecTrigName) != std::string::npos)
            //        {
            //            indexElecTrigger = i;
            //        }
            //        else if(triggerNames[i].find(metmhtTrigName) != std::string::npos)
            //        {
            //            indexMETMHTTrigger = i;
            //        }
            //    }
            //}
            //if(indexMuTrigger != -1 && indexElecTrigger != -1)
            //{
            //    // Check if the event passes the trigger, and float check that we are looking at the right trigger
            //    //if(triggerNames[indexMuTrigger].find(muTrigName) != std::string::npos && passTrigger[indexMuTrigger])
            //    //    passMuTrigger = true;
            //    if(triggerNames[indexElecTrigger].find(elecTrigName) != std::string::npos && passTrigger[indexElecTrigger])
            //        passElecTrigger = true;
            //    if(triggerNames[indexMETMHTTrigger].find(metmhtTrigName) != std::string::npos && passTrigger[indexMETMHTTrigger])
            //        passMETMHTTrigger = true;
            //}
            //else
            //{
            //    std::cout << "Could not find trigger in the list of trigger names" << std::endl;
            //}

            bool passSearchTrigger = false, passHighHtTrigger = false, passPhotonTrigger = false;
            for(int it = 0; it < triggerNames.size(); ++it)
            {
                if( triggerNames[it].find("HLT_PFMET170_NoiseCleaned_v")             != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET170_JetIdCleaned_v")             != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET170_HBHECleaned_v")              != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET100_PFMHT100_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET110_PFMHT110_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET120_PFMHT120_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET130_PFMHT130_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET140_PFMHT140_IDTight_v")         != std::string::npos || 
                    triggerNames[it].find("HLT_PFMET150_PFMHT150_IDTight_v")         != std::string::npos ||
                    triggerNames[it].find("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v") != std::string::npos ||
                    triggerNames[it].find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != std::string::npos ||
                    triggerNames[it].find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") != std::string::npos
                    )
                {
                    if( passTrigger[it] ) 
                    {
                        passSearchTrigger = true;
                    }
                }

                if( triggerNames[it].find("HLT_PFHT750_4JetPt50_v") != std::string::npos ||
                    triggerNames[it].find("HLT_PFHT800_v")          != std::string::npos ||
                    triggerNames[it].find("HLT_PFHT900_v")          != std::string::npos ||
                    triggerNames[it].find("HLT_PFJet450_v")         != std::string::npos
                    )
                {
                    if( passTrigger[it] ) 
                    {
                        passHighHtTrigger = true;
                    }                    
                }

                if( triggerNames[it].find("HLT_Photon175_v")                != std::string::npos ||
                    triggerNames[it].find("HLT_Photon75_v")                 != std::string::npos ||
                    triggerNames[it].find("HLT_Photon90_CaloIdL_PFHT500_v") != std::string::npos ||
                    triggerNames[it].find("HLT_Photon90_v")                 != std::string::npos
                    )
                {
                    if( passTrigger[it] ) 
                    {
                        passPhotonTrigger = true;
                    }                    
                }

                if( triggerNames[it].find("HLT_IsoMu24_v")              != std::string::npos ||
                    triggerNames[it].find("HLT_IsoTkMu24_v")            != std::string::npos ||
                    triggerNames[it].find("HLT_Mu50_v")                 != std::string::npos ||
                    triggerNames[it].find("HLT_Mu55_v")                 != std::string::npos
                    )
                {
                    if( passTrigger[it] ) 
                    {
                        passMuTrigger = true;
                    }                    
                }
            }

            tr.registerDerivedVar("passMuTrigger",     passMuTrigger);
            tr.registerDerivedVar("passElecTrigger",   passElecTrigger);
            tr.registerDerivedVar("passMETMHTTrigger", passMETMHTTrigger);
            tr.registerDerivedVar("passSearchTrigger", passSearchTrigger);
            tr.registerDerivedVar("passHighHtTrigger", passHighHtTrigger);
            tr.registerDerivedVar("passPhotonTrigger", passPhotonTrigger);
        }

        void triggerInfoMC(NTupleReader& tr)
        {
            const float& met                            = tr.getVar<float>("met");
            const float& ht                             = tr.getVar<float>("HT");
            const std::vector<TLorentzVector>& cutMuVec  = tr.getVec<TLorentzVector>("muonsLVec");

	    // MC trigger efficiencies
	    float triggerEff = GetTriggerEffWeight(met,ht);
	    float triggerEffStatUncUp = GetTriggerEffStatUncHi(met,ht);
	    float triggerEffSystUncUp = GetTriggerEffSystUncHi(met,ht);
	    float triggerEffUncUp     = TMath::Sqrt(triggerEffStatUncUp*triggerEffStatUncUp + triggerEffSystUncUp*triggerEffSystUncUp);
	    float triggerEffStatUncDown = GetTriggerEffStatUncLo(met,ht);
	    float triggerEffSystUncDown = GetTriggerEffSystUncLo(met,ht);
	    float triggerEffUncDown     = TMath::Sqrt(triggerEffStatUncDown*triggerEffStatUncDown + triggerEffSystUncDown*triggerEffSystUncDown);

            //Calculate muon trigger weights
            float muTrigWgt = 0.0;
            if(cutMuVec.size() >= 2 && cutMuVec[0].Pt() > 50 && cutMuVec[1].Pt() > 50)
            {
                float muEff1 = GetMuonTriggerEff(cutMuVec[0].Eta());
                float muEff2 = GetMuonTriggerEff(cutMuVec[1].Eta());

                muTrigWgt = 1 - (1 - muEff1)*(1 - muEff2);
            }
            else if(cutMuVec.size() >= 1 && cutMuVec[0].Pt() > 40)
            {
                //For events with only 1 muon (emu events in particular or events with a subleading muon below 45 GeV) just use the single muon eff
                muTrigWgt = GetMuonTriggerEff(cutMuVec[0].Eta());
            }

	    tr.registerDerivedVar("TriggerEffMC",triggerEff);
	    tr.registerDerivedVar("TriggerEffUpMC",triggerEff+triggerEffUncUp);
	    tr.registerDerivedVar("TriggerEffDownMC",triggerEff-triggerEffUncDown);

            tr.registerDerivedVar("muTrigWgt", muTrigWgt);
        }

    public:
	TriggerInfo(bool miniTuple = false, bool noMC = false)
	{
	    indexMuTrigger = -1;
	    indexElecTrigger = -1;
            indexMETMHTTrigger = -1;
            miniTuple_ = miniTuple;
            noMC_ = noMC;
	}

        void setIsMC(bool isMC)
        {
            if(isMC)
            {
                miniTuple_ = true;
                noMC_ = false;
            }
            else
            {
                miniTuple_ = false;
                noMC_ = true;
            }
        }

	void operator()(NTupleReader& tr)
	{
	    if(!miniTuple_) triggerInfo(tr);
            if(!noMC_)      triggerInfoMC(tr);
	}

    };

    class SystematicPrep
    {
    private:

        void systematicPrep(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& jetsLVec  = tr.getVec<TLorentzVector>("jetsLVec");
            const std::vector<float>& recoJetsJecUnc    = tr.getVec<float>("recoJetsJecUnc");

            const std::vector<float>& metMagUp   = tr.getVec<float>("metMagUp");
            const std::vector<float>& metMagDown = tr.getVec<float>("metMagDown");
            const std::vector<float>& metPhiUp   = tr.getVec<float>("metPhiUp");
            const std::vector<float>& metPhiDown = tr.getVec<float>("metPhiDown");

            const float& met    = tr.getVar<float>("met");
            const float& metphi = tr.getVar<float>("metphi");

            std::vector<TLorentzVector> *jetLVecUp = new std::vector<TLorentzVector>;
            std::vector<TLorentzVector> *jetLVecDn = new std::vector<TLorentzVector>;

            std::vector<float> *dPtMet = new std::vector<float>;
            std::vector<float> *dPhiMet = new std::vector<float>;

            float metUp = 0.0, metDn = 99990.0;

            for(int iMet = 0; iMet < metMagUp.size(); ++iMet)
            {
                metUp = std::max(metUp, metMagUp[iMet]);
                metDn = std::min(metDn, metMagDown[iMet]);
                
                dPtMet->push_back((metMagUp[iMet] - met)/met);
                dPtMet->push_back((metMagDown[iMet] - met)/met);
                dPhiMet->push_back(TVector2::Phi_mpi_pi(metPhiUp[iMet] - metphi));
                dPhiMet->push_back(TVector2::Phi_mpi_pi(metPhiDown[iMet] - metphi));
            }

            for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
            {
                jetLVecUp->push_back(jetsLVec[iJet] * (1 + recoJetsJecUnc[iJet]));
                jetLVecDn->push_back(jetsLVec[iJet] * (1 - recoJetsJecUnc[iJet]));
            }

            tr.registerDerivedVar("metMEUUp", metUp);
            tr.registerDerivedVar("metMEUDn", metDn);

            tr.registerDerivedVec("dPtMet", dPtMet);
            tr.registerDerivedVec("dPhiMet", dPhiMet);

            tr.registerDerivedVec("jetLVecUp", jetLVecUp);
            tr.registerDerivedVec("jetLVecDn", jetLVecDn);
        }

    public:
	SystematicPrep()
	{

	}

	void operator()(NTupleReader& tr)
	{
	    systematicPrep(tr);
	}

    };

    class SystematicCalc
    {
    private:
        SearchBins sbins;

        void systematicCalc(NTupleReader& tr)
        {
            const int& cntCSVSJEUUp = tr.getVar<int>("cntCSVSZinvJEUUp");
            const int& nTopCandSortedCntJEUUp = tr.getVar<int>("nTopCandSortedCntZinvJEUUp");
            const float& MT2JEUUp = tr.getVar<float>("best_had_brJet_MT2ZinvJEUUp");

            const int& cntCSVSJEUDn = tr.getVar<int>("cntCSVSZinvJEUDn");
            const int& nTopCandSortedCntJEUDn = tr.getVar<int>("nTopCandSortedCntZinvJEUDn");
            const float& MT2JEUDn = tr.getVar<float>("best_had_brJet_MT2ZinvJEUDn");

            const float& cleanMet = tr.getVar<float>("cleanMetPt");

            const int& cntCSVSMEUUp = tr.getVar<int>("cntCSVSZinvMEUUp");
            const int& nTopCandSortedCntMEUUp = tr.getVar<int>("nTopCandSortedCntZinvMEUUp");
            const float& MT2MEUUp = tr.getVar<float>("best_had_brJet_MT2ZinvMEUUp");
            const float& cleanMetMEUUp = tr.getVar<float>("metMEUUp");

            const int& cntCSVSMEUDn = tr.getVar<int>("cntCSVSZinvMEUDn");
            const int& nTopCandSortedCntMEUDn = tr.getVar<int>("nTopCandSortedCntZinvMEUDn");
            const float& MT2MEUDn = tr.getVar<float>("best_had_brJet_MT2ZinvMEUDn");
            const float& cleanMetMEUDn = tr.getVar<float>("metMEUDn");

            const float& HTUp           = tr.getVar<float>("HTZinvJEUUp");
            const float& HTDn           = tr.getVar<float>("HTZinvJEUDn");
            const float& HTMEUUp           = tr.getVar<float>("HTZinvMEUUp");
            const float& HTMEUDn           = tr.getVar<float>("HTZinvMEUDn");

            int nSearchBinJEUUp = sbins.find_Binning_Index(cntCSVSJEUUp, nTopCandSortedCntJEUUp, MT2JEUUp, cleanMet, HTUp);
            int nSearchBinJEUDn = sbins.find_Binning_Index(cntCSVSJEUDn, nTopCandSortedCntJEUDn, MT2JEUDn, cleanMet, HTDn);

            int nSearchBinMEUUp = sbins.find_Binning_Index(cntCSVSMEUUp, nTopCandSortedCntMEUUp, MT2MEUUp, cleanMetMEUUp, HTMEUUp);
            int nSearchBinMEUDn = sbins.find_Binning_Index(cntCSVSMEUDn, nTopCandSortedCntMEUDn, MT2MEUDn, cleanMetMEUDn, HTMEUDn);
            
            tr.registerDerivedVar("nSearchBinJEUUp", nSearchBinJEUUp);
            tr.registerDerivedVar("nSearchBinJEUDn", nSearchBinJEUDn);

            tr.registerDerivedVar("nSearchBinMEUUp", nSearchBinMEUUp);
            tr.registerDerivedVar("nSearchBinMEUDn", nSearchBinMEUDn);
        }

    public:
        SystematicCalc(std::string sb_era) : sbins(sb_era)
	{
            
	}

	void operator()(NTupleReader& tr)
	{
	    systematicCalc(tr);
	}

    };

    class PrepareTopCRSelection
    {
    private:

        int JECSys = 0;


        bool passNoiseEventFilterFunc(const NTupleReader* const tr, bool isfastsim = false)
        {
            // According to https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendationsICHEP16#Filters_to_be_applied,
            // "Do not apply filters to signal monte carlo (fastsim)"
            if( isfastsim ) return true;

            try
            {
                bool passDataSpec = true;
                if( tr->getVar<unsigned int>("run") >= 100000 ){ // hack to know if it's data or MC...
                    int goodVerticesFilter = tr->getVar<int>("goodVerticesFilter");
                    // new filters
                    const int & globalTightHalo2016Filter = tr->getVar<int>("globalTightHalo2016Filter");
                    bool passglobalTightHalo2016Filter = (&globalTightHalo2016Filter) != nullptr? tr->getVar<int>("globalTightHalo2016Filter") !=0 : true;

                    int eeBadScFilter = tr->getVar<int>("eeBadScFilter");

                    passDataSpec = goodVerticesFilter && eeBadScFilter && passglobalTightHalo2016Filter;
                }

                unsigned int hbheNoiseFilter = isfastsim? 1:tr->getVar<unsigned int>("HBHENoiseFilter");
                unsigned int hbheIsoNoiseFilter = isfastsim? 1:tr->getVar<unsigned int>("HBHEIsoNoiseFilter");
                int ecalTPFilter = tr->getVar<int>("EcalDeadCellTriggerPrimitiveFilter");

#ifdef TUPLE_OLD
                int jetIDFilter = isfastsim? 1:tr->getVar<int>("looseJetID");
#endif
#ifdef TUPLE_NEW
                int jetIDFilter = isfastsim? 1:tr->getVar<int>("AK4looseJetID");
#endif
                // new filters
                const unsigned int & BadPFMuonFilter = tr->getVar<unsigned int>("BadPFMuonFilter");
                bool passBadPFMuonFilter = (&BadPFMuonFilter) != nullptr? tr->getVar<unsigned int>("BadPFMuonFilter") !=0 : true;

                const unsigned int & BadChargedCandidateFilter = tr->getVar<unsigned int>("BadChargedCandidateFilter");
                bool passBadChargedCandidateFilter = (&BadChargedCandidateFilter) != nullptr? tr->getVar<unsigned int>("BadChargedCandidateFilter") !=0 : true;

                bool passMETratioFilter = tr->getVar<float>("calomet")!=0 ? tr->getVar<float>("met")/tr->getVar<float>("calomet") < 5 : true;

                //std::cout << (passDataSpec ? " TRUE" : "FALSE") << " " << (hbheNoiseFilter ? " TRUE" : "FALSE") << " " << (hbheIsoNoiseFilter ? " TRUE" : "FALSE") << " "
                //          << (ecalTPFilter ? " TRUE" : "FALSE") << " " << (jetIDFilter ? " TRUE" : "FALSE") << " " << (passBadPFMuonFilter ? " TRUE" : "FALSE") << " "
                //          << (passBadChargedCandidateFilter ? " TRUE" : "FALSE") << " " << (passMETratioFilter ? " TRUE" : "FALSE") <<std::endl;

                return passDataSpec && hbheNoiseFilter && hbheIsoNoiseFilter && ecalTPFilter && jetIDFilter && passBadPFMuonFilter && passBadChargedCandidateFilter && passMETratioFilter;
            }
            catch (std::string var)
            {
                if(tr->isFirstEvent()) 
                {
                    printf("NTupleReader::getTupleObj(const std::string var):  Variable not found: \"%s\"!!!\n", var.c_str());
                    printf("Running with PHYS14 Config\n");
                }
            }
            return true;
        }


        void prepareTopCRSelection(NTupleReader& tr)
        {
            //We need to to handle JEC systematics here.

            bool handleSys = true;

            const std::vector<TLorentzVector>& jetsLVecTemp = tr.getVec<TLorentzVector>("jetsLVec");
            const std::vector<float>& recoJetsJecUnc = tr.getVec<float>("recoJetsJecUnc");

            if(jetsLVecTemp.size() != recoJetsJecUnc.size()) handleSys = false; //If this is data, we can't do anything

            std::vector<TLorentzVector> jetsLVec;
            for(int ijet=0; ijet<jetsLVecTemp.size(); ++ijet)
            {
                if(handleSys){jetsLVec.push_back( jetsLVecTemp[ijet] * (1 + (JECSys * recoJetsJecUnc[ijet])));}
                else {jetsLVec.push_back( jetsLVecTemp[ijet] );}
            }            

//            const std::vector<TLorentzVector>& jetsLVec  = tr.getVec<TLorentzVector>("jetsLVec");
#ifdef TUPLE_OLD  
            //const std::vector<float>& recoJetsBtag      = tr.getVec<float>("recoJetsBtag_0");
#endif
#ifdef TUPLE_NEW  
            const std::vector<float>& recoJetsBtag      = tr.getVec<float>("recoJetsCSVv2");
#endif
	    const float& stored_weight = tr.getVar<float>("stored_weight");

            int cntCSVS = AnaFunctions::countCSVS(jetsLVec, recoJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);

            float bestCSV = -999;
            for(int i = 0; i < recoJetsBtag.size(); i++){
                bestCSV = (bestCSV > recoJetsBtag[i] ? bestCSV : recoJetsBtag[i]);
            }
            tr.registerDerivedVar("bestCSV",bestCSV);

            const float& metphi = tr.getVar<float>("metphi");

#ifdef TUPLE_NEW
            const std::vector<TLorentzVector>& gammaLVec = tr.getVec<TLorentzVector>("photonLVec");
#endif
#ifdef TUPLE_OLD
            const std::vector<TLorentzVector>& gammaLVec = tr.getVec<TLorentzVector>("gammaLVec");
#endif
            const std::vector<int>& tightPhotonID = tr.getVec<int>("tightPhotonID");

            std::vector<TLorentzVector> *tightPhotons = new std::vector<TLorentzVector>();

            //std::cout << "Comparing the length of the tightPhotonID and gammaLVec vectors: " << tightPhotonID.size() << " " << gammaLVec.size() << std::endl;

            int sizeP = (tightPhotonID.size() < gammaLVec.size() ? tightPhotonID.size() : gammaLVec.size());

            for(int i = 0; i < sizeP; ++i)
            {
                if(tightPhotonID[i])
                {
                    tightPhotons->push_back(gammaLVec[i]);
                }
            }

            tr.registerDerivedVec("tightPhotons", tightPhotons);
            tr.registerDerivedVar("passPhoton200", (tightPhotons->size() > 0) && ((*tightPhotons)[0].Pt() > 200));

            const std::vector<TLorentzVector>& muonsLVec    = tr.getVec<TLorentzVector>("muonsLVec");
            //const std::vector<float>& muonsRelIso          = tr.getVec<float>("muonsRelIso");
            const std::vector<float>& muonsMiniIso         = tr.getVec<float>("muonsMiniIso");
            const std::vector<float>& muonsMTlep           = tr.getVec<float>("muonsMtw");
            std::string muonsFlagIDLabel = "muonsFlagMedium";
            const std::vector<int> & muonsFlagIDVec = muonsFlagIDLabel.empty()? std::vector<int>(muonsMiniIso.size(), 1):tr.getVec<int>(muonsFlagIDLabel.c_str());

            std::vector<TLorentzVector>* cutMuVec = new std::vector<TLorentzVector>();
            std::vector<float> *cutMuMTlepVec = new std::vector<float>();
            for(int i = 0; i < muonsLVec.size(); ++i)
            {
                if(AnaFunctions::passMuon( muonsLVec[i], muonsMiniIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                {
                    cutMuVec->push_back(muonsLVec[i]);
                    cutMuMTlepVec->push_back(muonsMTlep[i]);
                }
            }

            tr.registerDerivedVec("cutMuVec", cutMuVec);
            tr.registerDerivedVec("cutMuMTlepVec", cutMuMTlepVec);

            const std::vector<TLorentzVector, std::allocator<TLorentzVector> > elesLVec = tr.getVec<TLorentzVector>("elesLVec");
            const std::vector<float>& elesMiniIso          = tr.getVec<float>("elesMiniIso");
            const std::vector<float>& elesCharge           = tr.getVec<float>("elesCharge");
            const std::vector<unsigned int>& elesisEB       = tr.getVec<unsigned int>("elesisEB",true);
            const std::vector<float>&  elesMTlep           = tr.getVec<float>("elesMtw");
            std::string elesFlagIDLabel = "elesFlagVeto";
            const std::vector<int> & elesFlagIDVec = elesFlagIDLabel.empty()? std::vector<int>(elesMiniIso.size(), 1):tr.getVec<int>(elesFlagIDLabel.c_str());

            //electron selection
            std::vector<TLorentzVector>* cutElecVec = new std::vector<TLorentzVector>();
            std::vector<float> *cutElecMTlepVec = new std::vector<float>();
            for(int i = 0; i < elesLVec.size(); ++i)
            {
                if(AnaFunctions::passElectron(elesLVec[i], elesMiniIso[i], -1, elesisEB[i], elesFlagIDVec[i], AnaConsts::elesMiniIsoArr))
                {
                    cutElecVec->push_back(elesLVec[i]);
                    cutElecMTlepVec->push_back(elesMTlep[i]);
                }
            }

            tr.registerDerivedVec("cutElecVec", cutElecVec);
            tr.registerDerivedVec("cutElecMTlepVec", cutElecMTlepVec);

            //// Calculate number of leptons
            int nMuons = AnaFunctions::countMuons(muonsLVec, muonsMiniIso, muonsMTlep, muonsFlagIDVec, AnaConsts::muonsMiniIsoArr);
            const AnaConsts::IsoAccRec muonsMiniIsoArr20GeV = {   -1,       2.4,      20,     -1,       0.2,     -1  };
            int nMuons_20GeV = AnaFunctions::countMuons(muonsLVec, muonsMiniIso, muonsMTlep, muonsFlagIDVec, muonsMiniIsoArr20GeV);
            const AnaConsts::IsoAccRec muonsMiniIsoArr30GeV = {   -1,       2.4,      30,     -1,       0.2,     -1  };
            int nMuons_30GeV = AnaFunctions::countMuons(muonsLVec, muonsMiniIso, muonsMTlep, muonsFlagIDVec, muonsMiniIsoArr30GeV);
            const AnaConsts::IsoAccRec muonsMiniIsoArr40GeV = {   -1,       2.4,      40,     -1,       0.2,     -1  };
            int nMuons_40GeV = AnaFunctions::countMuons(muonsLVec, muonsMiniIso, muonsMTlep, muonsFlagIDVec, muonsMiniIsoArr40GeV);
            const AnaConsts::IsoAccRec muonsMiniIsoArr50GeV = {   -1,       2.4,      50,     -1,       0.2,     -1  };
            int nMuons_50GeV = AnaFunctions::countMuons(muonsLVec, muonsMiniIso, muonsMTlep, muonsFlagIDVec, muonsMiniIsoArr50GeV);
            int nElectrons = AnaFunctions::countElectrons(elesLVec, elesMiniIso, elesMTlep, elesisEB,  elesFlagIDVec, AnaConsts::elesMiniIsoArr);
            const AnaConsts::ElecIsoAccRec elesMiniIsoArr20 = {   -1,       2.5,      20,     -1,     0.10,     0.10,     -1  };
            int nElectrons20 = AnaFunctions::countElectrons(elesLVec, elesMiniIso, elesMTlep, elesisEB, elesFlagIDVec, elesMiniIsoArr20);
            const AnaConsts::ElecIsoAccRec elesMiniIsoArr30 = {   -1,       2.5,      30,     -1,     0.10,     0.10,     -1  };
            int nElectrons30 = AnaFunctions::countElectrons(elesLVec, elesMiniIso, elesMTlep, elesisEB, elesFlagIDVec, elesMiniIsoArr30);
            int nIsoTrks; 
            if( tr.checkBranch("loose_isoTrksLVec") )
            {
                nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<float>("loose_isoTrks_iso"), tr.getVec<float>("loose_isoTrks_mtw"), tr.getVec<int>("loose_isoTrks_pdgId"));
            }
            else
            {
                nIsoTrks = 0;
            }
            //
            //// Pass lepton veto?
            bool passMuonVeto = (nMuons == AnaConsts::nMuonsSel);
            bool passEleVeto = (nElectrons == AnaConsts::nElectronsSel);
            bool passIsoTrkVeto = (nIsoTrks == AnaConsts::nIsoTrksSel);

            float Mmumu = -999.9;
            bool passfloatMuon = false;
            if(cutMuVec->size() >= 2)
            {
                Mmumu = ((*cutMuVec)[0] + (*cutMuVec)[1]).M();
                passfloatMuon = (*cutMuVec)[0].Pt() > 30 && (*cutMuVec)[1].Pt() > 20 && Mmumu > 81 && Mmumu < 101;
            }


            // Calculate deltaPhi
            std::vector<float> * dPhiVec = new std::vector<float>();
            (*dPhiVec) = AnaFunctions::calcDPhi(jetsLVec, metphi, 3, AnaConsts::dphiArr);

            // Pass deltaPhi?
            bool passdPhis = (dPhiVec->size() >= 3) && ((*dPhiVec)[0] >= AnaConsts::dPhi0_CUT && (*dPhiVec)[1] >= AnaConsts::dPhi1_CUT && (*dPhiVec)[2] >= AnaConsts::dPhi2_CUT);

            // calculate number of jets 
            int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);

            //calculate HT
            float HT = AnaFunctions::calcHT(jetsLVec, AnaConsts::pt30Eta24Arr);

	    // Process the generator weight
	    float genWeight = 1.;
	    // Never apply this weight for data! In the old ntuple version <=3 this is "-1", in the newer ones it is "0"
	    if(stored_weight < 0) genWeight = -1.;

            //std::cout << genWeight << std::endl;
	    tr.registerDerivedVar("genWeight", genWeight);

            tr.registerDerivedVar("cntCSVS", cntCSVS);

            tr.registerDerivedVar("passSingleLep50", nMuons_50GeV == 1);
            tr.registerDerivedVar("passSingleLep20", nMuons_20GeV + nElectrons20 == 1);
            tr.registerDerivedVar("passLep20", nMuons_20GeV + nElectrons20 >= 1);
            tr.registerDerivedVar("passSingleMu30", nMuons_30GeV == 1);
            tr.registerDerivedVar("passSingleMu40", nMuons_40GeV == 1);
            tr.registerDerivedVar("passSingleLep30", nMuons_30GeV + nElectrons30 == 1);
            tr.registerDerivedVar("passfloatLep", passfloatMuon);

            tr.registerDerivedVar("passLeptVetoNoMu", passEleVeto && passIsoTrkVeto);
            tr.registerDerivedVar("passLeptVeto", passMuonVeto && passEleVeto && passIsoTrkVeto);

            tr.registerDerivedVec("dPhiVec", dPhiVec);
            tr.registerDerivedVar("passdPhis", passdPhis);

            tr.registerDerivedVar("HT", HT);
            tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);
            
            tr.registerDerivedVar("passNoiseEventFilter", passNoiseEventFilterFunc(&tr));
        }

    public:
        void operator()(NTupleReader& tr)
        {
            prepareTopCRSelection(tr);
        }

        PrepareTopCRSelection(int JECcorr = 0){
            //Let's get ready for JEC systematics, we will set a variable here, and the Jet collection is loaded we will look to see what to do.
            if(JECcorr == -1){ JECSys = -1; }
            else if(JECcorr == 1){ JECSys = 1;}
            else{ JECSys = 0; } // If an invalid JECcorr value is pass, we will default to the regular behavior.
        }

    };
/*
    class AliasStealthVars
    {
    private:
        template<typename I, typename O> void addAliasVecTypeChange(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const std::vector<I>& inVec = tr.getVec<I>(name);
            std::vector<O>* outVec      = new std::vector<O>(inVec.size());
            for(int i = 0; i < inVec.size(); i++)
            {
                (*outVec)[i] = static_cast<O>(inVec[i]);
            }
            tr.registerDerivedVec(alias, outVec);
        }

        template<typename I, typename O> void addAliasTypeChange(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const I& inVar = tr.getVar<I>(name);
            O outVar = static_cast<O>(inVar);
            tr.registerDerivedVar(alias, outVar);
        }

        void inECALBarrel(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const std::vector<TLorentzVector>& vec = tr.getVec<TLorentzVector>(name);
            std::vector<unsigned int>* inBarrel    = new std::vector<unsigned int>(vec.size());
            for(int i = 0; i < vec.size(); i++)
            {
                TLorentzVector object = vec[i];
                double eta = object.Eta();
                if(fabs(eta < 1.479))
                {
                    (*inBarrel)[i] = static_cast<unsigned int>(1);
                }
                else
                {                  
                    (*inBarrel)[i] = static_cast<unsigned int>(0); 
                }
            }
            tr.registerDerivedVec(alias, inBarrel);
        }

        void genIDHack(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const std::vector<int>& vec = tr.getVec<int>(name);
            std::vector<int>* IdVec = new std::vector<int>(vec.size());
            for(int i = 0; i < vec.size(); i++)
            {
                (*IdVec)[i] = static_cast<int>(i);
            }
            tr.registerDerivedVec(alias, IdVec);
        }

        void add2Vec(NTupleReader& tr, const std::string& name1, const std::string& name2, const std::string& alias)
        {
            const std::vector<double>& vec1 = tr.getVec<double>(name1);
            const std::vector<double>& vec2 = tr.getVec<double>(name2);
            std::vector<double>* sumVec = new std::vector<double>(vec1.size());
            for(int i = 0; i < vec1.size(); i++)
            {
                (*sumVec)[i] = vec1[i] + vec2[i];
            }
            tr.registerDerivedVec(alias, sumVec);
        }

        void aliasVars(NTupleReader& tr)
        {
            inECALBarrel(tr,"Electrons","elesisEB");
            genIDHack(tr,"GenParticles_ParentIdx","genDecayIdxVec");
            add2Vec(tr,"Jets_neutralEmEnergyFraction","Jets_neutralHadronEnergyFraction","recoJetsneutralEnergyFraction");
            addAliasVecTypeChange<bool,int>(tr,"Muons_tightID","muonsFlagMedium");
            addAliasVecTypeChange<bool,int>(tr,"Electrons_tightID","elesFlagVeto");
            addAliasVecTypeChange<bool,int>(tr,"Photons_fullID","tightPhotonID");
            addAliasVecTypeChange<int,double>(tr,"Jets_chargedHadronMultiplicity","ChargedHadronMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_muonMultiplicity","MuonMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_neutralMultiplicity","NeutralHadronMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_photonMultiplicity","PhotonMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_electronMultiplicity","ElectronMultiplicity");
            addAliasTypeChange<bool,int>(tr,"JetID","looseJetID");
        }
        
    public:
        void addAllAlias(NTupleReader& tr)
        {
            tr.addAlias("Jets","jetsLVec");
            tr.addAlias("MET","met");
            tr.addAlias("Jets_bDiscriminatorCSV","recoJetsBtag_0");
            tr.addAlias("Weight","stored_weight");
            tr.addAlias("METPhi","metphi");
            tr.addAlias("GenParticles","genDecayLVec");
            tr.addAlias("puWeight","_PUweightFactor");
            tr.addAlias("Muons","muonsLVec");
            tr.addAlias("Muons_MiniIso","muonsMiniIso");
            tr.addAlias("Muons_MTW","muonsMtw");
            tr.addAlias("Electrons","elesLVec");
            tr.addAlias("Electrons_MiniIso","elesMiniIso");
            tr.addAlias("Electrons_charge","elesCharge");
            tr.addAlias("Electrons_MTW","elesMtw");
            tr.addAlias("Jets_qgLikelihood","qgLikelihood");
            tr.addAlias("Jets_ptD","qgPtD");
            tr.addAlias("Jets_axismajor","qgAxis1");
            tr.addAlias("Jets_axisminor","qgAxis2");
            tr.addAlias("Jets_chargedHadronEnergyFraction","recoJetschargedHadronEnergyFraction");
            tr.addAlias("Jets_chargedEmEnergyFraction","recoJetschargedEmEnergyFraction");
            tr.addAlias("Jets_neutralEmEnergyFraction","recoJetsneutralEmEnergyFraction");
            tr.addAlias("Jets_muonEnergyFraction","recoJetsmuonEnergyFraction");
            tr.addAlias("Jets_hfHadronEnergyFraction","recoJetsHFHadronEnergyFraction");
            tr.addAlias("Jets_hfEMEnergyFraction","recoJetsHFEMEnergyFraction");
            tr.addAlias("Jets_photonEnergyFraction","PhotonEnergyFraction");
            tr.addAlias("Jets_electronEnergyFraction","ElectronEnergyFraction");
            tr.addAlias("Jets_multiplicity","qgMult");
            tr.addAlias("Photons","gammaLVec");
            tr.addAlias("RunNum","run");
            //tr.addAlias("BadPFMuonFilter","BadPFMuonFilter");
            //tr.addAlias("BadChargedCandidateFilter","BadChargedCandidateFilter");
            tr.addAlias("CaloMET","calomet");
            tr.addAlias("GenParticles_PdgId","genDecayPdgIdVec");
            tr.addAlias("GenParticles_ParentIdx","genDecayMomIdxVec");
            tr.addAlias("TriggerPass","PassTrigger");
            tr.addAlias("Jets_partonFlavor","recoJetsFlavor"); //Could be Jets_hadronFlavor should ask Kevin
            tr.addAlias("NVtx","vtxSize");
        }

        void operator()(NTupleReader& tr)
        {
            aliasVars(tr);
        }
    };
*/
    class AliasStealthVars
    {
    private:
        template<typename I, typename O> void addAliasVecTypeChange(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const std::vector<I>& inVec = tr.getVec<I>(name);
            std::vector<O>* outVec      = new std::vector<O>(inVec.size());
            for(int i = 0; i < inVec.size(); i++)
            {
                (*outVec)[i] = static_cast<O>(inVec[i]);
            }
            tr.registerDerivedVec(alias, outVec);
        }

        template<typename I, typename O> void addAliasTypeChange(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const I& inVar = tr.getVar<I>(name);
            O outVar = static_cast<O>(inVar);
            tr.registerDerivedVar(alias, outVar);
        }

        void inECALBarrel(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const std::vector<TLorentzVector>& vec = tr.getVec<TLorentzVector>(name);
            std::vector<unsigned int>* inBarrel    = new std::vector<unsigned int>(vec.size());
            for(int i = 0; i < vec.size(); i++)
            {
                TLorentzVector object = vec[i];
                double eta = object.Eta();
                if(fabs(eta < 1.479))
                {
                    (*inBarrel)[i] = static_cast<unsigned int>(1);
                }
                else
                {                  
                    (*inBarrel)[i] = static_cast<unsigned int>(0); 
                }
            }
            tr.registerDerivedVec(alias, inBarrel);
        }

        void genIDHack(NTupleReader& tr, const std::string& name, const std::string& alias)
        {
            const std::vector<int>& vec = tr.getVec<int>(name);
            std::vector<int>* IdVec = new std::vector<int>(vec.size());
            for(int i = 0; i < vec.size(); i++)
            {
                (*IdVec)[i] = static_cast<int>(i);
            }
            tr.registerDerivedVec(alias, IdVec);
        }

        void add2Vec(NTupleReader& tr, const std::string& name1, const std::string& name2, const std::string& alias)
        {
            const std::vector<double>& vec1 = tr.getVec<double>(name1);
            const std::vector<double>& vec2 = tr.getVec<double>(name2);
            std::vector<double>* sumVec = new std::vector<double>(vec1.size());
            for(int i = 0; i < vec1.size(); i++)
            {
                (*sumVec)[i] = vec1[i] + vec2[i];
            }
            tr.registerDerivedVec(alias, sumVec);
        }

        void aliasVars(NTupleReader& tr)
        {
            inECALBarrel(tr,"Electrons","elesisEB");
            genIDHack(tr,"GenParticles_ParentIdx","genDecayIdxVec");
            add2Vec(tr,"Jets_neutralEmEnergyFraction","Jets_neutralHadronEnergyFraction","recoJetsneutralEnergyFraction");
            addAliasVecTypeChange<bool,int>(tr,"Muons_tightID","muonsFlagMedium");
            addAliasVecTypeChange<bool,int>(tr,"Electrons_tightID","elesFlagVeto");
            addAliasVecTypeChange<bool,int>(tr,"Photons_fullID","tightPhotonID");
            addAliasVecTypeChange<int,double>(tr,"Jets_chargedHadronMultiplicity","ChargedHadronMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_muonMultiplicity","MuonMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_neutralMultiplicity","NeutralHadronMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_photonMultiplicity","PhotonMultiplicity");
            addAliasVecTypeChange<int,double>(tr,"Jets_electronMultiplicity","ElectronMultiplicity");
            addAliasTypeChange<bool,int>(tr,"JetID","looseJetID");
        }
        
    public:
        void addAllAlias(NTupleReader& tr)
        {
            tr.addAlias("Jets","jetsLVec");
            tr.addAlias("MET","met");
            //tr.addAlias("Jets_bDiscriminatorCSV","recoJetsBtag_0");
            tr.addAlias("Weight","stored_weight");
            tr.addAlias("METPhi","metphi");
            tr.addAlias("GenParticles","genDecayLVec");
            tr.addAlias("puWeight","_PUweightFactor");
            tr.addAlias("Muons","muonsLVec");
            tr.addAlias("Muons_MiniIso","muonsMiniIso");
            tr.addAlias("Muons_MTW","muonsMtw");
            tr.addAlias("Electrons","elesLVec");
            tr.addAlias("Electrons_MiniIso","elesMiniIso");
            tr.addAlias("Electrons_charge","elesCharge");
            tr.addAlias("Electrons_MTW","elesMtw");
            tr.addAlias("Jets_qgLikelihood","qgLikelihood");
            tr.addAlias("Jets_ptD","qgPtD");
            tr.addAlias("Jets_axismajor","qgAxis1");
            tr.addAlias("Jets_axisminor","qgAxis2");
            tr.addAlias("Jets_chargedHadronEnergyFraction","recoJetschargedHadronEnergyFraction");
            tr.addAlias("Jets_chargedEmEnergyFraction","recoJetschargedEmEnergyFraction");
            tr.addAlias("Jets_neutralEmEnergyFraction","recoJetsneutralEmEnergyFraction");
            tr.addAlias("Jets_muonEnergyFraction","recoJetsmuonEnergyFraction");
            tr.addAlias("Jets_hfHadronEnergyFraction","recoJetsHFHadronEnergyFraction");
            tr.addAlias("Jets_hfEMEnergyFraction","recoJetsHFEMEnergyFraction");
            tr.addAlias("Jets_photonEnergyFraction","PhotonEnergyFraction");
            tr.addAlias("Jets_electronEnergyFraction","ElectronEnergyFraction");
            tr.addAlias("Jets_multiplicity","qgMult");
            tr.addAlias("Photons","gammaLVec");
            tr.addAlias("RunNum","run");
            //tr.addAlias("BadPFMuonFilter","BadPFMuonFilter");
            //tr.addAlias("BadChargedCandidateFilter","BadChargedCandidateFilter");
            tr.addAlias("CaloMET","calomet");
            tr.addAlias("GenParticles_PdgId","genDecayPdgIdVec");
            tr.addAlias("GenParticles_ParentIdx","genDecayMomIdxVec");
            tr.addAlias("TriggerPass","PassTrigger");
            tr.addAlias("Jets_partonFlavor","recoJetsFlavor"); //Could be Jets_hadronFlavor should ask Kevin
            tr.addAlias("NVtx","vtxSize");
        }

        void operator()(NTupleReader& tr)
        {
            aliasVars(tr);
        }
    };

    class PrepareTopVars
    {
    private:

        int indexMuTrigger, indexElecTrigger, indexHTMHTTrigger, indexMuHTTrigger;
        int JECSys = 0;
        std::shared_ptr<TopTagger> ttMVA;
        TopCat topMatcher_;
        std::shared_ptr<TFile> WMassCorFile;
        std::shared_ptr<TF1> puppisd_corrGEN;
        std::shared_ptr<TF1> puppisd_corrRECO_cen;
        std::shared_ptr<TF1> puppisd_corrRECO_for;
        
        std::mt19937 generator;
        std::uniform_int_distribution<int> distribution;

        void prepareTopVars(NTupleReader& tr)
        {
            //We need to to handle JEC systematics here.

            bool handleSys = true;

            const std::vector<TLorentzVector>& jetsLVecTemp = tr.getVec<TLorentzVector>("jetsLVec");
            const std::vector<float>& recoJetsJecUnc = tr.getVec<float>("recoJetsJecUnc");

            if(jetsLVecTemp.size() != recoJetsJecUnc.size()) handleSys = false; //If this is data, we can't do anything

            //std::cout << "handleSys is " << handleSys << std::endl;
            //std::cout << "JECSys is " << JECSys << std::endl;

            std::vector<TLorentzVector> jetsLVec;
            for(int ijet=0; ijet<jetsLVecTemp.size(); ++ijet)
            {
                if(handleSys){
                    jetsLVec.push_back( jetsLVecTemp[ijet] * (1 + (JECSys * recoJetsJecUnc[ijet])));
                    //std::cout << "Reco Jets uncertainty " << recoJetsJecUnc[ijet] << std::endl;
                }else{jetsLVec.push_back( jetsLVecTemp[ijet] );}
            }            

            for(int i = 0; i < jetsLVecTemp.size(); i++){
                //std::cout << "Jet pT before JEC unc: " << jetsLVecTemp[i].Pt() << ", Jet pT after JEC unc: " << jetsLVec[i].Pt() << std::endl;
            }

#ifdef TUPLE_OLD
            //const std::vector<float>& recoJetsBtag      = tr.getVec<float>("recoJetsBtag_0");
#endif
#ifdef TUPLE_NEW
            const std::vector<float>& recoJetsBtag      = tr.getVec<float>("recoJetsCSVv2");
#endif
            const std::vector<float>& qgLikelihood      = tr.getVec<float>("qgLikelihood");
            
            //AK8 variables 
            //const std::vector<float>& puppitau1    = tr.getVec<float>("puppitau1");
            //const std::vector<float>& puppitau2    = tr.getVec<float>("puppitau2");
            //const std::vector<float>& puppitau3    = tr.getVec<float>("puppitau3");
            //const std::vector<float>& puppisoftDropMass = tr.getVec<float>("puppisoftDropMass");
            //const std::vector<TLorentzVector>& puppiJetsLVec  = tr.getVec<TLorentzVector>("puppiJetsLVec");
            //const std::vector<TLorentzVector>& puppiSubJetsLVec  = tr.getVec<TLorentzVector>("puppiSubJetsLVec");
            //const std::vector<float>& puppiSubJetsBdisc = tr.getVec<float>("puppiSubJetsBdisc");
            //const std::vector<float>& puppiSubJetstotalMult = tr.getVec<float>("puppiSubJetstotalMult");
            //const std::vector<float>& puppiSubJetsptD = tr.getVec<float>("puppiSubJetsptD");
            //const std::vector<float>& puppiSubJetsaxis1 = tr.getVec<float>("puppiSubJetsaxis1");
            //const std::vector<float>& puppiSubJetsaxis2 = tr.getVec<float>("puppiSubJetsaxis2");
                        

            //Helper function to turn int vectors into float vectors
            auto convertTofloatandRegister = [](NTupleReader& tr, const std::string& name)
            {
                const std::vector<int>& intVec = tr.getVec<int>(name);
                std::vector<float>* floatVec = new std::vector<float>(intVec.begin(), intVec.end());
                tr.registerDerivedVec(name+"ConvertedTofloat3", floatVec);
                return floatVec;
            };
            
            //New Tagger starts here
            ttUtility::ConstAK4Inputs<float> *myConstAK4Inputs = nullptr;
            //ttUtility::ConstAK8Inputs<float> *myConstAK8Inputs = nullptr;
            std::vector<TLorentzVector> *genTops;
            std::vector<std::vector<const TLorentzVector*>> hadGenTopDaughters;
            std::vector<Constituent> constituentsMVA;
            if(tr.checkBranch("genDecayLVec"))
            {
                const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
                const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");
                const std::vector<int>& genDecayIdxVec          = tr.getVec<int>("genDecayIdxVec");
                const std::vector<int>& genDecayMomIdxVec       = tr.getVec<int>("genDecayMomIdxVec");

                //prep input object (constituent) vector
                genTops = new std::vector<TLorentzVector>(ttUtility::GetHadTopLVec(genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec));
                for(const auto& top : *genTops)
                {
                    hadGenTopDaughters.push_back(ttUtility::GetTopdauLVec(top, genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec));
                }

                myConstAK4Inputs = new ttUtility::ConstAK4Inputs<float>(jetsLVec, recoJetsBtag, qgLikelihood, *genTops, hadGenTopDaughters);

                //myConstAK8Inputs = new ttUtility::ConstAK8Inputs(puppiJetsLVec, puppitau1, puppitau2, puppitau3, puppisoftDropMass, puppiSubJetsLVec, puppiSubJetsBdisc, puppiSubJetstotalMult, puppiSubJetsptD, puppiSubJetsaxis1, puppiSubJetsaxis2, *genTops, hadGenTopDaughters);
            }
            else
            {
                //no gen info is avaliable
                genTops = new std::vector<TLorentzVector>();
                
                myConstAK4Inputs = new ttUtility::ConstAK4Inputs<float>(jetsLVec, recoJetsBtag, qgLikelihood);
                
                //myConstAK8Inputs = new ttUtility::ConstAK8Inputs(puppiJetsLVec, puppitau1, puppitau2, puppitau3, puppisoftDropMass, puppiSubJetsLVec, puppiSubJetsBdisc, puppiSubJetstotalMult, puppiSubJetsptD, puppiSubJetsaxis1, puppiSubJetsaxis2);
                
            }
                
            myConstAK4Inputs->addSupplamentalVector("qgLikelihood",                         tr.getVec<float>("qgLikelihood"));
            myConstAK4Inputs->addSupplamentalVector("qgPtD",                                tr.getVec<float>("qgPtD"));
            myConstAK4Inputs->addSupplamentalVector("qgAxis1",                              tr.getVec<float>("qgAxis1"));
            myConstAK4Inputs->addSupplamentalVector("qgAxis2",                              tr.getVec<float>("qgAxis2"));
            myConstAK4Inputs->addSupplamentalVector("recoJetschargedHadronEnergyFraction",  tr.getVec<float>("recoJetschargedHadronEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("recoJetschargedEmEnergyFraction",      tr.getVec<float>("recoJetschargedEmEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsneutralEmEnergyFraction",      tr.getVec<float>("recoJetsneutralEmEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsmuonEnergyFraction",           tr.getVec<float>("recoJetsmuonEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsHFHadronEnergyFraction",       tr.getVec<float>("recoJetsHFHadronEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsHFEMEnergyFraction",           tr.getVec<float>("recoJetsHFEMEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsneutralEnergyFraction",        tr.getVec<float>("recoJetsneutralEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("PhotonEnergyFraction",                 tr.getVec<float>("PhotonEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("ElectronEnergyFraction",               tr.getVec<float>("ElectronEnergyFraction"));
            myConstAK4Inputs->addSupplamentalVector("ChargedHadronMultiplicity",            tr.getVec<float>("ChargedHadronMultiplicity"));
            myConstAK4Inputs->addSupplamentalVector("NeutralHadronMultiplicity",            tr.getVec<float>("NeutralHadronMultiplicity"));
            myConstAK4Inputs->addSupplamentalVector("PhotonMultiplicity",                   tr.getVec<float>("PhotonMultiplicity"));
            myConstAK4Inputs->addSupplamentalVector("ElectronMultiplicity",                 tr.getVec<float>("ElectronMultiplicity"));
            myConstAK4Inputs->addSupplamentalVector("MuonMultiplicity",                     tr.getVec<float>("MuonMultiplicity"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVb",                             tr.getVec<float>("DeepCSVb"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVc",                             tr.getVec<float>("DeepCSVc"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVl",                             tr.getVec<float>("DeepCSVl"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVbb",                            tr.getVec<float>("DeepCSVbb"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVcc",                            tr.getVec<float>("DeepCSVcc"));
//            myConstAK4Inputs->addSupplamentalVector("DeepFlavorb",                          tr.getVec<float>("DeepFlavorb"));
//            myConstAK4Inputs->addSupplamentalVector("DeepFlavorbb",                         tr.getVec<float>("DeepFlavorbb"));
//            myConstAK4Inputs->addSupplamentalVector("DeepFlavorlepb",                       tr.getVec<float>("DeepFlavorlepb"));
//            myConstAK4Inputs->addSupplamentalVector("DeepFlavorc",                          tr.getVec<float>("DeepFlavorc"));
//            myConstAK4Inputs->addSupplamentalVector("DeepFlavoruds",                        tr.getVec<float>("DeepFlavoruds"));
//            myConstAK4Inputs->addSupplamentalVector("DeepFlavorg",                          tr.getVec<float>("DeepFlavorg"));
            myConstAK4Inputs->addSupplamentalVector("CvsL",                                 tr.getVec<float>("CversusL"));
//            myConstAK4Inputs->addSupplamentalVector("CvsB",                                 tr.getVec<float>("CvsB"));
            //myConstAK4Inputs->addSupplamentalVector("CombinedSvtx",                         tr.getVec<float>("CombinedSvtx"));
            //myConstAK4Inputs->addSupplamentalVector("JetProba",                             tr.getVec<float>("JetProba_0"));
            //myConstAK4Inputs->addSupplamentalVector("JetBprob",                             tr.getVec<float>("JetBprob"));
            //myConstAK4Inputs->addSupplamentalVector("recoJetsBtag",                         tr.getVec<float>("recoJetsBtag_0"));
            //myConstAK4Inputs->addSupplamentalVector("recoJetsCharge",                       tr.getVec<float>("recoJetsCharge_0"));
            myConstAK4Inputs->addSupplamentalVector("qgMult",                               *convertTofloatandRegister(tr, "qgMult"));
//            myConstAK4Inputs->addSupplamentalVector("qgMult",                               tr.getVec<float>("qgMult"));

            //myConstAK8Inputs.setWMassCorrHistos(puppisd_corrGEN, puppisd_corrRECO_cen, puppisd_corrRECO_for);

            //run new tagger
            //New MVA resolved Tagger starts here
/*
            int run = tr.getVar<int>("run");
            int lumi = tr.getVar<int>("lumi");
            long event = tr.getVar<long>("event");
 
            std::cout << std::endl << run << ":" << lumi << ":" << event << " (run:lumi:event) " << std::endl;
*/

            constituentsMVA = ttUtility::packageConstituents(*myConstAK4Inputs);//, *myConstAK8Inputs);
            //run tagger
            ttMVA->runTagger(constituentsMVA);

            delete myConstAK4Inputs;
            //delete myConstAK8Inputs;

            const TopTaggerResults& ttrMVA = ttMVA->getResults();
            //std::cout << "Top Candidates: " << ttrMVA.getTopCandidates().size() << std::endl;
            const auto& candidateTops = ttrMVA.getTopCandidates();
/*            for(int i = 0; i < candidateTops.size(); i++){
                auto& top = candidateTops[i];

                auto& j1 = top.getConstituents()[0]->p();
                auto& j2 = top.getConstituents()[1]->p();
                auto& j3 = top.getConstituents()[2]->p();

                std::cout << "Top Candidate #" << i << " disc: " << top.getDiscriminator() << ", eta: " << top.p().Eta() << ", phi: " << top.p().Phi() << ", pT: " << top.p().Pt() << ", cand_m: " << top.p().M()
                          << ", j12_m: " << (j1+j2).M() << ", j13_m: " << (j1+j3).M() << ", j23_m: " << (j2+j3).M() << std::endl;
                for(int j = 0; j < candidateTops[i].getConstituents().size(); j++){
                    auto& jet = top.getConstituents()[j];
                    auto qgAxis1 = jet->getExtraVar("qgAxis1");
                    auto qgAxis2 = jet->getExtraVar("qgAxis2");
                    auto qgMult = jet->getExtraVar("qgMult");
                    auto qgPtD = jet->getExtraVar("qgPtD");
                    auto ChargedHadronMultiplicity = jet->getExtraVar("ChargedHadronMultiplicity");
                    auto ElectronEnergyFraction = jet->getExtraVar("ElectronEnergyFraction");
                    auto ElectronMultiplicity = jet->getExtraVar("ElectronMultiplicity");
                    auto MuonMultiplicity = jet->getExtraVar("MuonMultiplicity");
                    auto NeutralHadronMultiplicity = jet->getExtraVar("NeutralHadronMultiplicity");
                    auto PhotonEnergyFraction = jet->getExtraVar("PhotonEnergyFraction");
                    auto PhotonMultiplicity = jet->getExtraVar("PhotonMultiplicity");
                    auto recoJetsHFEMEEnergyFraction = jet->getExtraVar("recoJetsHFEMEnergyFraction");
                    auto recoJetsHFHadronEnergyFraction = jet->getExtraVar("recoJetsHFHadronEnergyFraction");
                    auto recoJetschargedEmEnergyFraction = jet->getExtraVar("recoJetschargedEmEnergyFraction");
                    auto recoJetschargedHadronEnergyFraction = jet->getExtraVar("recoJetschargedHadronEnergyFraction");
                    auto recoJetsmuonEnergyFraction = jet->getExtraVar("recoJetsmuonEnergyFraction");
                    auto recoJetsneutralEmEnergyFraction = jet->getExtraVar("recoJetsneutralEmEnergyFraction");
                    auto recoJetsneutralEnergyFraction = jet->getExtraVar("recoJetsneutralEnergyFraction");
                    auto DeepCSVb = jet->getExtraVar("DeepCSVb");
                    auto DeepCSVbb = jet->getExtraVar("DeepCSVbb");
                    auto DeepCSVc = jet->getExtraVar("DeepCSVc");
                    auto DeepCSVcc = jet->getExtraVar("DeepCSVcc");
                    auto DeepCSVl = jet->getExtraVar("DeepCSVl");
                    std::cout << "Jet #" << j+1 << ") m:" << jet->p().M() << ", p: " << jet->p().P()
                              << ", qgAxis1: " << qgAxis1
                              << ", qgAxis2: " << qgAxis2
                              << ", qgMult: " << qgMult
                              << ", qgPtD: " << qgPtD
                              << ", ChargedHadronMultiplicity: " << ChargedHadronMultiplicity
                              << ", ElectronEnergyFraction: " << ElectronEnergyFraction
                              << ", ElectronMultiplicity: " << ElectronMultiplicity
                              << ", MuonMultiplicity: " << MuonMultiplicity
                              << ", NeutralHadronMultiplicity: " << NeutralHadronMultiplicity
                              << ", PhotonEnergyFraction: " << PhotonEnergyFraction
                              << ", PhotonMultiplicity: " << PhotonMultiplicity
                              << ", recoJetsHFEMEEnergyFraction: " << recoJetsHFEMEEnergyFraction
                              << ", recoJetsHFHadronEnergyFraction: " << recoJetsHFHadronEnergyFraction
                              << ", recoJetschargedEmEnergyFraction: " << recoJetschargedEmEnergyFraction
                              << ", recoJetschargedHadronEnergyFraction: " << recoJetschargedHadronEnergyFraction
                              << ", recoJetsmuonEnergyFraction: " << recoJetsmuonEnergyFraction
                              << ", recoJetsneutralEmEnergyFraction: " << recoJetsneutralEmEnergyFraction
                              << ", recoJetsneutralEnergyFraction: " << recoJetsneutralEnergyFraction
                              << ", DeepCSVb: " << DeepCSVb
                              << ", DeepCSVbb: " << DeepCSVbb
                              << ", DeepCSVc: " << DeepCSVc
                              << ", DeepCSVcc: " << DeepCSVcc
                              << ", DeepCSVl: " << DeepCSVl << std::endl;
                }
                std::cout << std::endl;
            }*/
            const auto& tops = ttrMVA.getTops();

            //get "best" top based upon on trijet mass 
            float bestTopMass = -9999.9;
            float bestTopEta = -9999.9;
            const TopObject* bestTopMassLV = nullptr;
            bool bestTopMassGenMatch = false;
            bool bestTopMassTopTag = false;
            
            float highestDisc = -9999.9;

            for(int iTop = 0; iTop < candidateTops.size(); ++iTop)
            {
                auto& top = candidateTops[iTop];

                highestDisc = (top.getDiscriminator() > highestDisc ? top.getDiscriminator() : highestDisc);

                if(fabs(top.p().M() - 173.5) < fabs(bestTopMass - 173.5) && top.getNConstituents() == 3)
                {
                    bestTopMass = top.p().M();
                    bestTopEta = top.p().Eta();
                    bestTopMassLV = &top;
                }
            }

            bestTopMassGenMatch = (bestTopMassLV)?(bestTopMassLV->getBestGenTopMatch(0.6) != nullptr):(false);
            for(const auto& topPtr : tops) 
            {
                if(topPtr == bestTopMassLV) 
                {
                    bestTopMassTopTag = true;
                    break;
                }
            }

            //get a random top candidate 
            const TopObject* randomTopCand = nullptr;
            bool randomTopCandTopTag = false;
            bool randomTopCandNotGenMatch = true;
            if(candidateTops.size())
            {
                randomTopCand = &candidateTops[distribution(generator) % candidateTops.size()];
                randomTopCandNotGenMatch = (randomTopCand->getBestGenTopMatch(0.6) == nullptr);
                for(const auto& topPtr : tops) 
                {
                    if(topPtr == randomTopCand) 
                    {
                        randomTopCandTopTag = true;
                        break;
                    }
                }
            }


            tr.registerDerivedVar("ttrMVA", &ttrMVA);

            //get one mu of 20 GeV pt
            tr.registerDerivedVar("nTops", static_cast<int>(tops.size()));

            tr.registerDerivedVec("genTops", genTops);

            tr.registerDerivedVar("highestDisc", highestDisc);

            tr.registerDerivedVar("bestTopMass", bestTopMass);
            tr.registerDerivedVar("bestTopEta", bestTopEta);
            tr.registerDerivedVar("bestTopMassLV", bestTopMassLV?(bestTopMassLV->p()):(TLorentzVector()));
            tr.registerDerivedVar("bestTopMassGenMatch", bestTopMassGenMatch);
            tr.registerDerivedVar("bestTopMassTopTag", bestTopMassTopTag);


            tr.registerDerivedVar("randomTopCand", randomTopCand?(randomTopCand->p()):(TLorentzVector()));
            tr.registerDerivedVar("randomTopCandNConst", randomTopCand?(randomTopCand->getNConstituents()):(-1));
            tr.registerDerivedVar("randomTopCandTopTag", randomTopCandTopTag);
            tr.registerDerivedVar("randomTopCandNotGenMatch", randomTopCandNotGenMatch);

            //std::cout << "Finished prepareTopVar" << std::endl;

        }


    public:
        PrepareTopVars(std::string taggerCfg = "TopTagger.cfg", int JECcorr = 0) : ttMVA(new TopTagger()), WMassCorFile(nullptr), puppisd_corrGEN(nullptr), puppisd_corrRECO_cen(nullptr), puppisd_corrRECO_for(nullptr), distribution(1,65000)
	{
            //Let's get ready for JEC systematics, we will set a variable here, and the Jet collection is loaded we will look to see what to do.
            if(JECcorr == -1){ JECSys = -1; }
            else if(JECcorr == 1){ JECSys = 1;}
            else{ JECSys = 0; } // If an invalid JECcorr value is pass, we will default to the regular behavior.

            ttMVA->setCfgFile(taggerCfg);

            indexMuTrigger = indexElecTrigger = indexHTMHTTrigger = indexMuHTTrigger = -1;

            std::string puppiCorr = "puppiCorr.root";
            WMassCorFile.reset(TFile::Open(puppiCorr.c_str(),"READ"));
            if (!WMassCorFile)
                std::cout << "W mass correction file not found w mass!!!!!!! " << puppiCorr <<" Will not correct W mass" << std::endl;
            else{
                puppisd_corrGEN     .reset((TF1*)WMassCorFile->Get("puppiJECcorr_gen"));
                puppisd_corrRECO_cen.reset((TF1*)WMassCorFile->Get("puppiJECcorr_reco_0eta1v3"));
                puppisd_corrRECO_for.reset((TF1*)WMassCorFile->Get("puppiJECcorr_reco_1v3eta2v5"));
            }

	}

        ~PrepareTopVars()
        {
            //if(tt) delete tt;
        }

	void operator()(NTupleReader& tr)
	{
            prepareTopVars(tr);
	}

    };


    class MetSmear
    {
    private:

        TRandom3 *trand;

        float calcB(float x1, float x2, float e1, float e2)
        {
            return (x2*x2*(e1-1) - x1*x1*(e2-1)) / (x1*x2*(x2-x1));
        }

        float calcC(float x1, float x2, float e1, float e2)
        {
            return (x2*(e1-1) - x1*(e2-1)) / (x1*x2*(x1-x2));
        }

        float calcQuad(float x, float x1, float x2, float e1, float e2)
        {
            return 1 + x*calcB(x1, x2, e1, e2) + x*x*calcC(x1, x2, e1, e2);
        }

        float logistical(float met, float A, float B, float C)
        {
            return 1 - A/(1 + exp(-(B*(met - C))));
        }
        
        void metSmear(NTupleReader& tr)
        {
            const float& met     = tr.getVar<float>("cleanMetPt");
            
            // Logistical smearing 
            float met_logi_1 = met * logistical(met, 0.15, 0.01, 300);
            float met_logi_2 = met * logistical(met, 0.20, 0.01, 400);
            float met_logi_3 = met * logistical(met, 0.25, 0.01, 500);
            float met_logi_4 = met * logistical(met, 0.20, 0.01, 400);
            float met_logi_5 = met * logistical(met, 0.20, 0.02, 400);
            float met_logi_6 = met * logistical(met, 0.20, 0.03, 400);
            float met_logi_7 = met * logistical(met, 0.20, 0.02, 300);
            float met_logi_8 = met * logistical(met, 0.20, 0.02, 400);
            float met_logi_9 = met * logistical(met, 0.20, 0.02, 500);

            // gaussian smearing 
            //float met_gaus_5  = trand->Gaus(met, 5);
            //float met_gaus_10 = trand->Gaus(met, 10);
            //float met_gaus_15 = trand->Gaus(met, 15);
            float met_gaus_20 = trand->Gaus(met, 20);
            //float met_gaus_25 = trand->Gaus(met, 25);
            float met_gaus_30 = trand->Gaus(met, 30);
            float met_gaus_40 = trand->Gaus(met, 40);
            float met_gaus_50 = trand->Gaus(met, 50);

            tr.registerDerivedVar("met_gaus_20", met_gaus_20);
            //tr.registerDerivedVar("met_gaus_25", met_gaus_25);
            tr.registerDerivedVar("met_gaus_30", met_gaus_30);
            tr.registerDerivedVar("met_gaus_40", met_gaus_40);
            tr.registerDerivedVar("met_gaus_50", met_gaus_50);

            tr.registerDerivedVar("met_logi_1", met_logi_1);
            tr.registerDerivedVar("met_logi_2", met_logi_2);
            tr.registerDerivedVar("met_logi_3", met_logi_3);
            tr.registerDerivedVar("met_logi_4", met_logi_4);
            tr.registerDerivedVar("met_logi_5", met_logi_5);
            tr.registerDerivedVar("met_logi_6", met_logi_6);
            tr.registerDerivedVar("met_logi_7", met_logi_7);
            tr.registerDerivedVar("met_logi_8", met_logi_8);
            tr.registerDerivedVar("met_logi_9", met_logi_9);
        }

        void mt2Smear(NTupleReader& tr)
        {
            const float& metphi       = tr.getVar<float>("cleanMetPhi");
            const float& met_logi_1   = tr.getVar<float>("met_logi_1");
            const float& met_gaus_30  = tr.getVar<float>("met_gaus_30");
            
            const std::vector<TLorentzVector>& jetsLVec_forTagger  = tr.getVec<TLorentzVector>("jetsLVec_forTaggerZinv");
            const std::vector<float>&     recoJetsBtag_forTagger  = tr.getVec<float>("recoJetsBtag_forTaggerZinv");

            //We choose 30 GeV gaussian smearing and logi_1 for the study

            // Form TLorentzVector of MET
            TLorentzVector metLVec_Logi;
            metLVec_Logi.SetPtEtaPhiM(met_logi_1, 0, metphi, 0);
            
            //type3Ptr->processEvent(jetsLVec_forTagger, recoJetsBtag_forTagger, metLVec_Logi);
            float MT2_Logi = 0.0;//type3Ptr->best_had_brJet_MT2;

            TLorentzVector metLVec_Gaus;
            metLVec_Gaus.SetPtEtaPhiM(met_gaus_30, 0, metphi, 0);
            
            //type3Ptr->processEvent(jetsLVec_forTagger, recoJetsBtag_forTagger, metLVec_Gaus); 
            float MT2_Gaus = 0.0;//type3Ptr->best_had_brJet_MT2;

            tr.registerDerivedVar("mt2_logi_1",  MT2_Logi);
            tr.registerDerivedVar("mt2_gaus_30", MT2_Gaus);
        }

    public:
        MetSmear()
        {
            trand = new TRandom3(452147);
        }

        //~MetSmear()
        // {
        //    if(trand) delete trand;
        // }

        void operator()(NTupleReader& tr)
        {
            metSmear(tr);
            mt2Smear(tr);
        }
    };

    class PrepareMiniTupleVars
    {
    private:

        static const int BIT_PASSLEPTVETO             = 0x00000001;
        static const int BIT_PASSMUONVETO             = 0x00000002;
        static const int BIT_PASSELEVETO              = 0x00000004;
        static const int BIT_PASSISOTRKVETO           = 0x00000008;
        static const int BIT_PASSNJETS                = 0x00000010;
        static const int BIT_PASSDPHIS                = 0x00000020;
        static const int BIT_PASSBJETS                = 0x00000040;
        static const int BIT_PASSMET                  = 0x00000080;
        static const int BIT_PASSMT2                  = 0x00000100;
        static const int BIT_PASSHT                   = 0x00000200;
        static const int BIT_PASSTAGGER               = 0x00000400;
        static const int BIT_PASSNOISEEVENTFILTER     = 0x00000800;
        static const int BIT_PASSBASELINE             = 0x00001000;
        static const int BIT_PASSBASELINENOTAGMT2     = 0x00002000;
        static const int BIT_PASSBASELINENOTAG        = 0x00004000;
        static const int BIT_PASSLEPTVETOZINV         = 0x00008000;
        static const int BIT_PASSMUONVETOZINV         = 0x00010000;
        static const int BIT_PASSELEVETOZINV          = 0x00020000;
        static const int BIT_PASSISOTRKVETOZINV       = 0x00040000;
        static const int BIT_PASSNJETSZINV            = 0x00080000;
        static const int BIT_PASSDPHISZINV            = 0x00100000;
        static const int BIT_PASSBJETSZINV            = 0x00200000;
        static const int BIT_PASSMETZINV              = 0x00400000;
        static const int BIT_PASSMT2ZINV              = 0x00800000;
        static const int BIT_PASSHTZINV               = 0x01000000;
        static const int BIT_PASSTAGGERZINV           = 0x02000000;
        static const int BIT_PASSNOISEEVENTFILTERZINV = 0x04000000;
        static const int BIT_PASSBASELINEZINV         = 0x08000000;
        static const int BIT_PASSBASELINENOTAGMT2ZINV = 0x10000000;
        static const int BIT_PASSBASELINENOTAGZINV    = 0x20000000;
        static const int BIT_PASSMUZINVSEL            = 0x40000000;
        static const int BIT_PASSELMUZINVSEL          = 0x80000000;

        bool pack_;

        void pack(NTupleReader& tr)
        {
            const bool& passLeptVeto =         tr.getVar<bool>("passLeptVeto");
            const bool& passMuonVeto =         tr.getVar<bool>("passMuonVeto");
            const bool& passEleVeto =          tr.getVar<bool>("passEleVeto");
            const bool& passIsoTrkVeto =       tr.getVar<bool>("passIsoTrkVeto");
            const bool& passnJets =            tr.getVar<bool>("passnJets");
            const bool& passdPhis =            tr.getVar<bool>("passdPhis");
            const bool& passBJets =            tr.getVar<bool>("passBJets");
            const bool& passMET =              tr.getVar<bool>("passMET");
            const bool& passMT2 =              tr.getVar<bool>("passMT2");
            const bool& passHT =               tr.getVar<bool>("passHT");
            const bool& passTagger =           tr.getVar<bool>("passTagger");
            const bool& passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilter");
            const bool& passBaseline =         tr.getVar<bool>("passBaseline");
            const bool& passBaselineNoTagMT2 = tr.getVar<bool>("passBaselineNoTagMT2");
            const bool& passBaselineNoTag =    tr.getVar<bool>("passBaselineNoTag");

            const bool& passLeptVetoZinv =         tr.getVar<bool>("passLeptVetoZinv");
            const bool& passMuonVetoZinv =         tr.getVar<bool>("passMuonVetoZinv");
            const bool& passEleVetoZinv =          tr.getVar<bool>("passEleVetoZinv");
            const bool& passIsoTrkVetoZinv =       tr.getVar<bool>("passIsoTrkVetoZinv");
            const bool& passnJetsZinv =            tr.getVar<bool>("passnJetsZinv");
            const bool& passdPhisZinv =            tr.getVar<bool>("passdPhisZinv");
            const bool& passBJetsZinv =            tr.getVar<bool>("passBJetsZinv");
            const bool& passMETZinv =              tr.getVar<bool>("passMETZinv");
            const bool& passMT2Zinv =              tr.getVar<bool>("passMT2Zinv");
            const bool& passHTZinv =               tr.getVar<bool>("passHTZinv");
            const bool& passTaggerZinv =           tr.getVar<bool>("passTaggerZinv");
            const bool& passNoiseEventFilterZinv = tr.getVar<bool>("passNoiseEventFilterZinv");
            const bool& passBaselineZinv =         tr.getVar<bool>("passBaselineZinv");
            const bool& passBaselineNoTagMT2Zinv = tr.getVar<bool>("passBaselineNoTagMT2Zinv");
            const bool& passBaselineNoTagZinv =    tr.getVar<bool>("passBaselineNoTagZinv");

            const bool& passMuZinvSel =            tr.getVar<bool>("passMuZinvSel");
            const bool& passElMuZinvSel =          tr.getVar<bool>("passElMuZinvSel");

            int cuts = 0;

            if(passLeptVeto)         cuts |= BIT_PASSLEPTVETO;
            if(passMuonVeto)         cuts |= BIT_PASSMUONVETO;
            if(passEleVeto)          cuts |= BIT_PASSELEVETO;
            if(passIsoTrkVeto)       cuts |= BIT_PASSISOTRKVETO;
            if(passnJets)            cuts |= BIT_PASSNJETS;
            if(passdPhis)            cuts |= BIT_PASSDPHIS;
            if(passBJets)            cuts |= BIT_PASSBJETS;
            if(passMET)              cuts |= BIT_PASSMET;
            if(passMT2)              cuts |= BIT_PASSMT2;
            if(passHT)               cuts |= BIT_PASSHT;
            if(passTagger)           cuts |= BIT_PASSTAGGER;
            if(passNoiseEventFilter) cuts |= BIT_PASSNOISEEVENTFILTER;
            if(passBaseline)         cuts |= BIT_PASSBASELINE;
            if(passBaselineNoTagMT2) cuts |= BIT_PASSBASELINENOTAGMT2;
            if(passBaselineNoTag)    cuts |= BIT_PASSBASELINENOTAG;

            if(passLeptVetoZinv)         cuts |= BIT_PASSLEPTVETOZINV;
            if(passMuonVetoZinv)         cuts |= BIT_PASSMUONVETOZINV;
            if(passEleVetoZinv)          cuts |= BIT_PASSELEVETOZINV;
            if(passIsoTrkVetoZinv)       cuts |= BIT_PASSISOTRKVETOZINV;
            if(passnJetsZinv)            cuts |= BIT_PASSNJETSZINV;
            if(passdPhisZinv)            cuts |= BIT_PASSDPHISZINV;
            if(passBJetsZinv)            cuts |= BIT_PASSBJETSZINV;
            if(passMETZinv)              cuts |= BIT_PASSMETZINV;
            if(passMT2Zinv)              cuts |= BIT_PASSMT2ZINV;
            if(passHTZinv)               cuts |= BIT_PASSHTZINV;
            if(passTaggerZinv)           cuts |= BIT_PASSTAGGERZINV;
            if(passNoiseEventFilterZinv) cuts |= BIT_PASSNOISEEVENTFILTERZINV;
            if(passBaselineZinv)         cuts |= BIT_PASSBASELINEZINV;
            if(passBaselineNoTagMT2Zinv) cuts |= BIT_PASSBASELINENOTAGMT2ZINV;
            if(passBaselineNoTagZinv)    cuts |= BIT_PASSBASELINENOTAGZINV;

            if(passMuZinvSel)            cuts |= BIT_PASSMUZINVSEL;
            if(passElMuZinvSel)          cuts |= BIT_PASSELMUZINVSEL;

            tr.registerDerivedVar("cuts", cuts);
        }

        void unpack(NTupleReader& tr)
        {
            const int& cuts = tr.getVar<int>("cuts");

            tr.registerDerivedVar("passLeptVeto",         static_cast<bool>(cuts & BIT_PASSLEPTVETO));
            tr.registerDerivedVar("passMuonVeto",         static_cast<bool>(cuts & BIT_PASSMUONVETO));
            tr.registerDerivedVar("passEleVeto",          static_cast<bool>(cuts & BIT_PASSELEVETO));
            tr.registerDerivedVar("passIsoTrkVeto",       static_cast<bool>(cuts & BIT_PASSISOTRKVETO));
            tr.registerDerivedVar("passnJets",            static_cast<bool>(cuts & BIT_PASSNJETS));
            tr.registerDerivedVar("passdPhis",            static_cast<bool>(cuts & BIT_PASSDPHIS));
            tr.registerDerivedVar("passBJets",            static_cast<bool>(cuts & BIT_PASSBJETS));
            tr.registerDerivedVar("passMET",              static_cast<bool>(cuts & BIT_PASSMET));
            tr.registerDerivedVar("passMT2",              static_cast<bool>(cuts & BIT_PASSMT2));
            tr.registerDerivedVar("passHT",               static_cast<bool>(cuts & BIT_PASSHT));
            tr.registerDerivedVar("passTagger",           static_cast<bool>(cuts & BIT_PASSTAGGER));
            tr.registerDerivedVar("passNoiseEventFilter", static_cast<bool>(cuts & BIT_PASSNOISEEVENTFILTER));
            tr.registerDerivedVar("passBaseline",         static_cast<bool>(cuts & BIT_PASSBASELINE));
            tr.registerDerivedVar("passBaselineNoTagMT2", static_cast<bool>(cuts & BIT_PASSBASELINENOTAGMT2));
            tr.registerDerivedVar("passBaselineNoTag",    static_cast<bool>(cuts & BIT_PASSBASELINENOTAG));

            tr.registerDerivedVar("passLeptVetoZinv",         static_cast<bool>(cuts & BIT_PASSLEPTVETOZINV));
            tr.registerDerivedVar("passMuonVetoZinv",         static_cast<bool>(cuts & BIT_PASSMUONVETOZINV));
            tr.registerDerivedVar("passEleVetoZinv",          static_cast<bool>(cuts & BIT_PASSELEVETOZINV));
            tr.registerDerivedVar("passIsoTrkVetoZinv",       static_cast<bool>(cuts & BIT_PASSISOTRKVETOZINV));
            tr.registerDerivedVar("passnJetsZinv",            static_cast<bool>(cuts & BIT_PASSNJETSZINV));
            tr.registerDerivedVar("passdPhisZinv",            static_cast<bool>(cuts & BIT_PASSDPHISZINV));
            tr.registerDerivedVar("passBJetsZinv",            static_cast<bool>(cuts & BIT_PASSBJETSZINV));
            tr.registerDerivedVar("passMETZinv",              static_cast<bool>(cuts & BIT_PASSMETZINV));
            tr.registerDerivedVar("passMT2Zinv",              static_cast<bool>(cuts & BIT_PASSMT2ZINV));
            tr.registerDerivedVar("passHTZinv",               static_cast<bool>(cuts & BIT_PASSHTZINV));
            tr.registerDerivedVar("passTaggerZinv",           static_cast<bool>(cuts & BIT_PASSTAGGERZINV));
            tr.registerDerivedVar("passNoiseEventFilterZinv", static_cast<bool>(cuts & BIT_PASSNOISEEVENTFILTERZINV));
            tr.registerDerivedVar("passBaselineZinv",         static_cast<bool>(cuts & BIT_PASSBASELINEZINV));
            tr.registerDerivedVar("passBaselineNoTagMT2Zinv", static_cast<bool>(cuts & BIT_PASSBASELINENOTAGMT2ZINV));
            tr.registerDerivedVar("passBaselineNoTagZinv",    static_cast<bool>(cuts & BIT_PASSBASELINENOTAGZINV));

            tr.registerDerivedVar("passMuZinvSel",   static_cast<bool>(cuts & BIT_PASSMUZINVSEL));
            tr.registerDerivedVar("passElMuZinvSel", static_cast<bool>(cuts & BIT_PASSELMUZINVSEL));
        }

    public:
        PrepareMiniTupleVars(bool pack)
        {
            pack_ = pack;
        }

        void operator()(NTupleReader& tr)
        {
            if(pack_) pack(tr);
            else      unpack(tr);
        }
    };

    class NJetAk8 {
     private:
        void generateNJetAk8(NTupleReader& tr) {
          const std::vector<TLorentzVector>& ak8JetsLVec  = tr.getVec<TLorentzVector>("ak8JetsLVec"); 
          const std::vector<TLorentzVector>& puppiJetsLVec  = tr.getVec<TLorentzVector>("puppiJetsLVec");
        }
     public:
     NJetAk8(){}
    ~NJetAk8(){}
             void operator()(NTupleReader& tr)
          {
            generateNJetAk8(tr);
          }
        };   


     class Taudiv {
      private:
          std::shared_ptr<TopTagger> ttPtr_mine;
          void generateTaudiv(NTupleReader& tr) {
            const std::vector<float>& tau1    = tr.getVec<float>("tau1");
            const std::vector<float>& tau2    = tr.getVec<float>("tau2");
            const std::vector<float>& tau3    = tr.getVec<float>("tau3");
            const std::vector<float>& puppitau1    = tr.getVec<float>("puppitau1");
            const std::vector<float>& puppitau2    = tr.getVec<float>("puppitau2");
            const std::vector<float>& puppitau3    = tr.getVec<float>("puppitau3");
            const std::vector<float>& softDropMass = tr.getVec<float>("softDropMass");
            const std::vector<float>& puppisoftDropMass = tr.getVec<float>("puppisoftDropMass");
            const std::vector<TLorentzVector>& jetsLVec     = tr.getVec<TLorentzVector>("jetsLVec");
            const std::vector<TLorentzVector>& ak8JetsLVec  = tr.getVec<TLorentzVector>("ak8JetsLVec");
            const std::vector<TLorentzVector>& puppiJetsLVec  = tr.getVec<TLorentzVector>("puppiJetsLVec");

            std::vector<TLorentzVector> *puppiLVecLoose_top = new std::vector<TLorentzVector>();
            std::vector<TLorentzVector> *puppiLVectight_top = new std::vector<TLorentzVector>();
            std::vector<TLorentzVector> *puppiLVecLoose_w = new std::vector<TLorentzVector>();
            std::vector<TLorentzVector> *puppiLVectight_w = new std::vector<TLorentzVector>();
            std::vector<float>* puppitau2Dtau1 = new std::vector<float>();
            std::vector<float>* puppitau3Dtau2 = new std::vector<float>();
            std::vector<float>* puppitau2Dtau1_SDM = new std::vector<float>();
            std::vector<float>* puppitau3Dtau2_SDM = new std::vector<float>();

            std::vector<TLorentzVector> *hadWLVec = new std::vector<TLorentzVector>();

            int monoJet=0;
            int diJet=0;
            int triJet=0;
            const TopTaggerResults& ttr =ttPtr_mine->getResults();
            std::vector<TopObject*> Ntop = ttr.getTops();
            for(int i=0; i<Ntop.size(); i++){
                if(Ntop[i]->getNConstituents() == 1) monoJet++;
                else if(Ntop[i]->getNConstituents() == 2) diJet++;
                else if(Ntop[i]->getNConstituents() == 3) triJet++;
            }
            
            const int& nJetsAk8 = ak8JetsLVec.size(); 
            const int& nJetsPuppi = puppiJetsLVec.size();
            tr.registerDerivedVar("nJetsAk8", nJetsAk8);
            tr.registerDerivedVar("nJetsPuppi", nJetsPuppi);
            tr.registerDerivedVar("typeMono",monoJet);
            tr.registerDerivedVar("typeDi",diJet);
            tr.registerDerivedVar("typeTri",triJet);
            
            if(puppitau2.size()!=0 && puppitau1.size()!=0 && puppitau2.size()==puppitau1.size()){
		for(int iJet = 0; iJet < nJetsPuppi; ++iJet){
		    puppitau2Dtau1->push_back(puppitau2[iJet]/(puppitau1[iJet]));
		}
            }
            else { 
		puppitau2Dtau1->push_back( -1);
	    }

            if(puppitau2.size()!=0 && puppitau3.size()!=0 && puppitau2.size()==puppitau3.size()){
		for(int iJet = 0; iJet < nJetsPuppi; ++iJet){
		    puppitau3Dtau2->push_back(puppitau3[iJet]/(puppitau2[iJet]));
		}
	    }
            else {
		puppitau3Dtau2->push_back( -1);
            }
	    tr.registerDerivedVec("puppitau2Dtau1", puppitau2Dtau1);
	    tr.registerDerivedVec("puppitau3Dtau2", puppitau3Dtau2);
            
	    ///WTagging
            for(int tau = 0; tau < (*puppitau2Dtau1).size(); ++tau){
               if (puppisoftDropMass[tau]>65 && puppisoftDropMass[tau]<100){
		   // push back tau variables after mass cut
		   puppitau2Dtau1_SDM->push_back(puppitau2Dtau1->at(tau));

		   if ((*puppitau2Dtau1)[tau] >= 0 && (*puppitau2Dtau1)[tau] < 0.6){ // loose
		       puppiLVecLoose_w->push_back(puppiJetsLVec[tau]);  
		       //std::cout <<"PT_puupi"<< (*puppiLVectight_w).size()  << std::endl;    // (*puppiLVectight_w)[0].Pt() 

		       if ((*puppitau2Dtau1)[tau] < 0.45){ // tight
			   puppiLVectight_w->push_back(puppiJetsLVec[tau]); 
		       }
		   }
               } 
	    }

            //Top 1%
	    for(int tau = 0; tau < (*puppitau3Dtau2).size(); ++tau){
		if (puppisoftDropMass[tau]>105 && puppisoftDropMass[tau]<210){
		    puppitau3Dtau2_SDM->push_back(puppitau3Dtau2->at(tau));
		    
		    if ((*puppitau3Dtau2)[tau] >= 0 && (*puppitau3Dtau2)[tau] < 0.65){
			puppiLVecLoose_top->push_back(puppiJetsLVec[tau]);

			if ((*puppitau3Dtau2)[tau] < 0.54){ 
			    puppiLVectight_top->push_back(puppiJetsLVec[tau]);
			}
		    }
		}
	    }
	    
	    tr.registerDerivedVec("puppiLVectight_top", puppiLVectight_top);
	    tr.registerDerivedVec("puppiLVecLoose_top", puppiLVecLoose_top);
	    tr.registerDerivedVec("puppiLVectight_w", puppiLVectight_w);
	    tr.registerDerivedVec("puppiLVecLoose_w", puppiLVecLoose_w);
	    tr.registerDerivedVec("puppitau2Dtau1_SDM", puppitau2Dtau1_SDM);
	    tr.registerDerivedVec("puppitau3Dtau2_SDM", puppitau3Dtau2_SDM);

	  }

        public:
          Taudiv(std::shared_ptr<TopTagger> ttPtr) { 
           ttPtr_mine = ttPtr;
          }
          ~Taudiv() {}
          void operator()(NTupleReader& tr)
          {
            generateTaudiv(tr);
          }       

        };

     class Ak8DrMatch {
     private:
	 void generateAk8DrMatch(NTupleReader& tr) {
	     const std::vector<TLorentzVector>& jetsLVec     = tr.getVec<TLorentzVector>("jetsLVec");
	     const std::vector<TLorentzVector>& ak8JetsLVec  = tr.getVec<TLorentzVector>("ak8JetsLVec");
	     const std::vector<TLorentzVector>& puppiLVectight_top = tr.getVec<TLorentzVector>("puppiLVectight_top");
	     const std::vector<TLorentzVector>& puppiLVecLoose_top = tr.getVec<TLorentzVector>("puppiLVecLoose_top");
	     const std::vector<TLorentzVector>& puppiLVectight_w = tr.getVec<TLorentzVector>("puppiLVectight_w");
	     const std::vector<TLorentzVector>& puppiLVecLoose_w = tr.getVec<TLorentzVector>("puppiLVecLoose_w");  
	     const std::vector<TLorentzVector>& puppiJetsLVec  = tr.getVec<TLorentzVector>("puppiJetsLVec");

	     int nJetsAK41_min = 0;
	     int nJetsAK41_med = 0; 
	     int nJetsAK41_lar = 0;
	     int nJetsPuppi_T1_min = 0;
	     int nJetsPuppi_T1_med = 0;
	     int nJetsPuppi_T1_lar = 0;
	     int nJetsPuppi_L1_min = 0;
	     int nJetsPuppi_L1_med = 0;
	     int nJetsPuppi_L1_lar = 0;
	     int nJetsAK42_min = 0;
	     int nJetsAK42_med = 0;
	     int nJetsAK42_lar = 0;
 
	     std::vector<float>* ak81dRMin = new std::vector<float>();
	     std::vector<float>* ak82dRMin = new std::vector<float>(); 
	     std::vector<float>* puppi_top_L_1dRMin = new std::vector<float>();
	     std::vector<float>* puppi_top_L_2dRMin = new std::vector<float>();
	     std::vector<float>* puppi_top_T_1dRMin = new std::vector<float>();
	     std::vector<float>* puppi_top_T_2dRMin = new std::vector<float>(); 
	     for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
	     {
		 if(ak8JetsLVec.size() >= 1) ak81dRMin->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], ak8JetsLVec[0]));
		 if(ak8JetsLVec.size() >= 2) ak82dRMin->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], ak8JetsLVec[1]));
	//	 if(puppiLVectight_top.size() >= 1) puppi_top_L_1dRMin-> push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVectight_top[0]));
	//	 if(puppiLVectight_top.size() >= 2) puppi_top_L_2dRMin ->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVectight_top[1]));
	//	 if(puppiLVecLoose_top.size() >= 1) puppi_top_T_1dRMin-> push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVecLoose_top[0]));
	//	 if(puppiLVecLoose_top.size() >= 2) puppi_top_T_2dRMin ->push_back( ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], puppiLVecLoose_top[1]));
		 std::sort( ak81dRMin->begin(),ak81dRMin->end() );
		 std::sort( ak82dRMin->begin(),ak82dRMin->end() );
	     }
	     tr.registerDerivedVec("ak81dRMin", ak81dRMin);
	     tr.registerDerivedVec("ak82dRMin", ak82dRMin);
	  //   tr.registerDerivedVec("puppi_top_L_1dRMin", puppi_top_L_1dRMin);
	   //  tr.registerDerivedVec("puppi_top_L_2dRMin", puppi_top_L_2dRMin);    
	    // tr.registerDerivedVec("puppi_top_T_1dRMin", puppi_top_T_1dRMin);
	     //tr.registerDerivedVec("puppi_top_T_2dRMin", puppi_top_T_2dRMin);     
 
	     for(int iJet1 = 0; iJet1 < ak81dRMin->size(); ++iJet1)
	     {
		 if(ak81dRMin->at(iJet1) <=0.2) nJetsAK41_min++;
		 if(ak81dRMin->at(iJet1) <=0.4 && ak81dRMin->at(iJet1) > 0.2) nJetsAK41_med++;
		 if(ak81dRMin->at(iJet1) <=0.8 && ak81dRMin->at(iJet1) > 0.4) nJetsAK41_lar++;
            
             } 
	     for(int iJet1 = 0; iJet1 < puppi_top_L_1dRMin->size(); ++iJet1)
	     {
		 if(puppi_top_L_1dRMin->at(iJet1) <=0.2) nJetsPuppi_L1_min++;
		 if(puppi_top_L_1dRMin->at(iJet1) <=0.4 && puppi_top_L_1dRMin->at(iJet1) > 0.2) nJetsPuppi_L1_med++;
		 if(puppi_top_L_1dRMin->at(iJet1) <=0.8 && puppi_top_L_1dRMin->at(iJet1) > 0.4) nJetsPuppi_L1_lar++;

             }
	     for(int iJet1 = 0; iJet1 < puppi_top_T_1dRMin->size(); ++iJet1)
	     {
		 if(puppi_top_T_1dRMin->at(iJet1) <=0.2) nJetsPuppi_T1_min++;
		 if(puppi_top_T_1dRMin->at(iJet1) <=0.4 && puppi_top_T_1dRMin->at(iJet1) > 0.2) nJetsPuppi_T1_med++;
		 if(puppi_top_T_1dRMin->at(iJet1) <=0.8 && puppi_top_T_1dRMin->at(iJet1) > 0.4) nJetsPuppi_T1_lar++;

             }
	     for(int iJet2 = 0; iJet2 < ak82dRMin->size(); ++iJet2)
	     {
		 if(ak82dRMin->at(iJet2) <=0.2) nJetsAK42_min++;
		 if(ak82dRMin->at(iJet2) <=0.4 && ak82dRMin->at(iJet2) > 0.2) nJetsAK42_med++;                         
		 if(ak82dRMin->at(iJet2) <=0.8 && ak82dRMin->at(iJet2) > 0.4) nJetsAK42_lar++;
   
	     }          
	     tr.registerDerivedVar("nJetsAK41_min",nJetsAK41_min);
	     tr.registerDerivedVar("nJetsAK41_med",nJetsAK41_med);
	     tr.registerDerivedVar("nJetsAK41_lar",nJetsAK41_lar);
	     tr.registerDerivedVar("nJetsPuppi_L1_min",nJetsPuppi_L1_min);
	     tr.registerDerivedVar("nJetsPuppi_L1_med",nJetsPuppi_L1_med);
	     tr.registerDerivedVar("nJetsPuppi_L1_lar",nJetsPuppi_L1_lar);
	     tr.registerDerivedVar("nJetsPuppi_T1_min",nJetsPuppi_T1_min);
	     tr.registerDerivedVar("nJetsPuppi_T1_med",nJetsPuppi_T1_med);
	     tr.registerDerivedVar("nJetsPuppi_T1_lar",nJetsPuppi_T1_lar);
	     tr.registerDerivedVar("nJetsAK42_min",nJetsAK42_min);
	     tr.registerDerivedVar("nJetsAK42_med",nJetsAK42_med);
	     tr.registerDerivedVar("nJetsAK42_lar",nJetsAK42_lar);


	     // Also start looking at subjet information
	     const std::vector<TLorentzVector>& puppiSubJetsLVec  = tr.getVec<TLorentzVector>("puppiSubJetsLVec");
	     const std::vector<float>& puppiSubJetsBdisc = tr.getVec<float>("puppiSubJetsBdisc");

	     // For each tagged top/W, find the corresponding subjets
	     /*
	     std::vector< std::vector<TLorentzVector> > W_subjets;
	     std::vector<float>* W_subjets_pt_reldiff = new std::vector<float>();
	     for( TLorentzVector myW : puppiLVectight_w)
	     {
		 std::vector<TLorentzVector> myW_subjets;
		 int i = 0;
		 for(TLorentzVector puppiSubJet : puppiSubJetsLVec)
		 {
		     float myDR = ROOT::Math::VectorUtil::DeltaR(myW, puppiSubJet);
		     if (myDR < 0.8)
		     {
			 myW_subjets.push_back(puppiSubJet);
		     }
		     ++i;
		 }
		 // If more than 2 matches, find the best combination of two subjets by checking diff in 4-vector
		 if (myW_subjets.size() > 2) {
		     float min_diff = 999999.;
		     int min_j=0, min_k=1;
		     for (int j=0 ; j<myW_subjets.size(); ++j)
		     {
			 for (int k=j+1; k<myW_subjets.size(); ++k)
			 {
			     TLorentzVector diff_LV = myW - myW_subjets[j] - myW_subjets[k];
			     float diff = abs(diff_LV.M());
			     if(diff < min_diff)
			     {
				 min_diff = diff;
				 min_j = j;
				 min_k = k;
			     }
			 }
		     }
		     std::vector<TLorentzVector> mynewW_subjets = {myW_subjets[min_j], myW_subjets[min_k]};
		     W_subjets.push_back(mynewW_subjets);
		     W_subjets_pt_reldiff->push_back( ((myW_subjets[min_j]+myW_subjets[min_k]).Pt()-myW.Pt())/myW.Pt());
		 } else {
		     W_subjets.push_back(myW_subjets);
		     W_subjets_pt_reldiff->push_back( ((myW_subjets[0]+myW_subjets[1]).Pt()-myW.Pt())/myW.Pt());
		 }
	     }
	     tr.registerDerivedVec("W_subjets_pt_reldiff", W_subjets_pt_reldiff);
             */
	     // For each tagged top/W, find the corresponding subjets
	     std::vector< std::vector< TLorentzVector> > top_subjets;
	     std::vector<float>* top_subjets_pt_reldiff = new std::vector<float>();
             /*
	     for( TLorentzVector mytop : puppiLVectight_top)
	     {
		 std::vector<TLorentzVector> mytop_subjets;
		 int i = 0;
		 for(TLorentzVector puppiSubJet : puppiSubJetsLVec)
		 {
		     float myDR = ROOT::Math::VectorUtil::DeltaR(mytop, puppiSubJet);
		     if (myDR < 0.8)
		     {
			 mytop_subjets.push_back(puppiSubJet);
		     }
		     ++i;
		 }
		 // If more than 2 matches, find the best combination of two subjets
		 if (mytop_subjets.size() > 2) {
		     float min_diff = 999999.;
		     int min_j=0, min_k=1;
		     for (int j=0 ; j<mytop_subjets.size(); ++j)
		     {
			 for (int k=j+1; k<mytop_subjets.size(); ++k)
			 {
			     TLorentzVector diff_LV = mytop - mytop_subjets[j] - mytop_subjets[k];
			     float diff = abs(diff_LV.M());
			     if(diff < min_diff)
			     {
				 min_diff = diff;
				 min_j = j;
				 min_k = k;
			     }
			 }
		     }
		     std::vector<TLorentzVector> mynewtop_subjets = {mytop_subjets[min_j], mytop_subjets[min_k]};
		     top_subjets.push_back(mynewtop_subjets);
		     top_subjets_pt_reldiff->push_back( ((mytop_subjets[min_j]+mytop_subjets[min_k]).Pt()-mytop.Pt())/mytop.Pt());
		 } else {
		     top_subjets.push_back(mytop_subjets);
		     top_subjets_pt_reldiff->push_back( ((mytop_subjets[0]+mytop_subjets[1]).Pt()-mytop.Pt())/mytop.Pt());
		 }
	     }
	     tr.registerDerivedVec("top_subjets_pt_reldiff", top_subjets_pt_reldiff);
             */
	     // Figure out gen matching..
	     std::vector<bool>* gentop_match = new std::vector<bool>(); // helpful to make plots of matched and unmatched number of tops
	     std::vector<float>* dR_top_gentop = new std::vector<float>(); 
	     std::vector<float>* dR_AK4_topsubjet_genmatched = new std::vector<float>(); 
	     std::vector<float>* dR_AK4_top_genmatched = new std::vector<float>(); 
	     std::vector<int>* top_N_AK4_matched_genmatched = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_matched_notgenmatched = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_notmatched_genmatched = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_notmatched_notgenmatched = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_matched_genmatched_0p6 = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_matched_notgenmatched_0p6 = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_notmatched_genmatched_0p6 = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_notmatched_notgenmatched_0p6 = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_matched_genmatchedother = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_matched_notgenmatchedother = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_matched_genmatchedother_0p6 = new std::vector<int>(); 
	     std::vector<int>* top_N_AK4_matched_notgenmatchedother_0p6 = new std::vector<int>(); 
	     if(tr.checkBranch("genDecayLVec") && tr.checkBranch("genDecayPdgIdVec"))
	     {
                 const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");
                 const std::vector<int>& genDecayIdxVec          = tr.getVec<int>("genDecayIdxVec");
                 const std::vector<int>& genDecayMomIdxVec       = tr.getVec<int>("genDecayMomIdxVec");
                 const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
		 // For each tagged top, find the matching gen particles

		 // These are the hadronically decaying top quarks in the event:
		 std::vector<TLorentzVector> hadtopLVec = genUtility::GetHadTopLVec(genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec);
		 std::vector< std::vector<TLorentzVector> > hadtopdauLVec;
		 for(TLorentzVector hadtop : hadtopLVec)
		 {
		     hadtopdauLVec.push_back(genUtility::GetTopdauLVec(hadtop, genDecayLVec, genDecayPdgIdVec, genDecayIdxVec, genDecayMomIdxVec));
		 }

		 // check all tagged tops
		 /*
		 for(unsigned int imytop=0; imytop<puppiLVectight_top.size(); ++imytop) 
		 {
		     TLorentzVector mytop = puppiLVectight_top[imytop];
		     //std::cout << "Mytop info: " << mytop.Pt() << " " << mytop.Eta() << " " << mytop.Phi() << std::endl;
		     // For now find the closest hadtop in deltaR
		     TLorentzVector temp_gentop_match_LV;
		     float min_DR = 99.;
		     int matched_hadtop_index = -1;
		     for(unsigned int myhadtop_i=0; myhadtop_i<hadtopLVec.size(); ++myhadtop_i)
		     {
			 TLorentzVector myhadtop = hadtopLVec[myhadtop_i];
			 float DR_top = ROOT::Math::VectorUtil::DeltaR(mytop, myhadtop);
			 if (DR_top < min_DR) 
			 {
			     temp_gentop_match_LV = myhadtop;
			     min_DR = DR_top;
			     matched_hadtop_index = myhadtop_i;
			 }
		     }
		     dR_top_gentop->push_back(min_DR);
		     // DR should be small for it to actually be a match
		     if(min_DR < 0.4)
		     {
			 //std::cout << "Mytop info: " << mytop.Pt() << " " << mytop.Eta() << " " << mytop.Phi() << std::endl;
			 gentop_match->push_back(true);
			 // Now find the gen daughters for this gentop
			 std::vector<TLorentzVector> gentopdauLVec = hadtopdauLVec[matched_hadtop_index];

			 // Now we have the tagged top (mytop), the gen had top (temp_gentop_match_LV), and the gen daughters (gentopdauLVec)
			 // ready for some matching FUN!

			 // Removing AK4 jets based on DR matching with subjets of tagged top
			 std::vector<TLorentzVector> mysubjets = top_subjets[imytop];
			 std::vector<int> ak4_removed;
			 if(mysubjets.size() != 2)
			     std::cout << "Attention: found " << mysubjets.size() << " subjets instead of 2" << std::endl;
			 //std::cout << "Subjet 0: " << mysubjets[0].Pt() << " " << mysubjets[0].Eta() << " " << mysubjets[0].Phi() << std::endl;
			 //std::cout << "Subjet 0: " << mysubjets[1].Pt() << " " << mysubjets[1].Eta() << " " << mysubjets[1].Phi() << std::endl;
			 
			 // some counters
			 int N_AK4_matched_genmatched = 0;
			 int N_AK4_matched_notgenmatched = 0;
			 int N_AK4_notmatched_genmatched = 0;
			 int N_AK4_notmatched_notgenmatched = 0;
			 // some counters
			 int N_AK4_matched_genmatched_0p6 = 0;
			 int N_AK4_matched_notgenmatched_0p6 = 0;
			 int N_AK4_notmatched_genmatched_0p6 = 0;
			 int N_AK4_notmatched_notgenmatched_0p6 = 0;
			 // matched to another gentop
			 int N_AK4_matched_genmatchedother = 0;
			 int N_AK4_matched_notgenmatchedother = 0;
			 int N_AK4_matched_genmatchedother_0p6 = 0;
			 int N_AK4_matched_notgenmatchedother_0p6 = 0;

			 for (unsigned int j=0; j<jetsLVec.size(); ++j)
			 {
			     float DR1 = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], mysubjets[0]);
			     float DR2 = DR1;
                             if(mysubjets.size()>1) {
                             float DR2 = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], mysubjets[1]);
                             }
			     //std::cout << "DR1, DR2: " << DR1 << " " << DR2 << std::endl;
			     // Check if it matches a gen daughter
			     bool genmatch = false;
			     for (TLorentzVector gendau : gentopdauLVec)
			     {
				 float DR_AK4_gen = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], gendau);
				 //std::cout << "gen DR " << DR_AK4_gen << std::endl;
				 if (DR_AK4_gen < 0.4)
				 {
				     // matches gendaughter
				     genmatch = true;
				     break;
				 }
			     }
			     if(genmatch){
				 dR_AK4_topsubjet_genmatched->push_back(std::min(DR1,DR2));
				 dR_AK4_top_genmatched->push_back(ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], mytop));
			     }
			     // should merge this with 'genmatch' finding...
			     bool genmatch_other = false;
			     for (unsigned int other=0; other<hadtopdauLVec.size() && !genmatch_other; ++other)
			     {
				 if(other == matched_hadtop_index)
				     continue;
				 for (TLorentzVector gendau : hadtopdauLVec[other])
				 {
				     float DR_AK4_gen = ROOT::Math::VectorUtil::DeltaR(jetsLVec[j], gendau);
				     //std::cout << "gen DR " << DR_AK4_gen << std::endl;
				     if (DR_AK4_gen < 0.4)
				     {
					 // matches gendaughter
					 genmatch_other = true;
					 break;
				     }
				 }
			     }

			     bool subjetmatch = false;
			     bool subjetmatch_0p6 = false;
			     if (DR1 < 0.4 || DR2 < 0.4)
			     {
				 //std::cout << "Found AK4 jet matching a subjet" << std::endl;
				 // found a match
				 subjetmatch = true;
				 ak4_removed.push_back(j);
			     }
			     if (DR1 < 0.6 || DR2 < 0.6)
			     {
				 //std::cout << "Found AK4 jet matching a subjet" << std::endl;
				 // found a match
				 subjetmatch_0p6 = true;
			     }


			     if(genmatch){
				 if(subjetmatch)
				     N_AK4_matched_genmatched++;
				 else
				     N_AK4_notmatched_genmatched++;
				 if(subjetmatch_0p6)
				     N_AK4_matched_genmatched_0p6++;
				 else
				     N_AK4_notmatched_genmatched_0p6++;
			     } else { // Not genmatched to any of the correct top daughters
				 if(subjetmatch)
				     N_AK4_matched_notgenmatched++;
				 else
				     N_AK4_notmatched_notgenmatched++;
				 if(subjetmatch_0p6)
				     N_AK4_matched_notgenmatched_0p6++;
				 else
				     N_AK4_notmatched_notgenmatched_0p6++;

				 if(genmatch_other){
				     if(subjetmatch)
					 N_AK4_matched_genmatchedother++;
				     if(subjetmatch_0p6)
					 N_AK4_matched_genmatchedother_0p6++;
				 } else {
				     if(subjetmatch)
					 N_AK4_matched_notgenmatchedother++;
				     if(subjetmatch_0p6)
					 N_AK4_matched_notgenmatchedother_0p6++;
				 }

			     }




			 }
			 top_N_AK4_matched_genmatched->push_back(N_AK4_matched_genmatched);
			 top_N_AK4_matched_notgenmatched->push_back(N_AK4_matched_notgenmatched);
			 top_N_AK4_notmatched_genmatched->push_back(N_AK4_notmatched_genmatched);
			 top_N_AK4_notmatched_notgenmatched->push_back(N_AK4_notmatched_notgenmatched);
			 top_N_AK4_matched_genmatched_0p6->push_back(N_AK4_matched_genmatched_0p6);
			 top_N_AK4_matched_notgenmatched_0p6->push_back(N_AK4_matched_notgenmatched_0p6);
			 top_N_AK4_notmatched_genmatched_0p6->push_back(N_AK4_notmatched_genmatched_0p6);
			 top_N_AK4_notmatched_notgenmatched_0p6->push_back(N_AK4_notmatched_notgenmatched_0p6);

			 top_N_AK4_matched_genmatchedother->push_back(N_AK4_matched_genmatchedother);
			 top_N_AK4_matched_notgenmatchedother->push_back(N_AK4_matched_notgenmatchedother);
			 top_N_AK4_matched_genmatchedother_0p6->push_back(N_AK4_matched_genmatchedother_0p6);
			 top_N_AK4_matched_notgenmatchedother_0p6->push_back(N_AK4_matched_notgenmatchedother_0p6);
			 
		     } else // No match
		     { 
			 gentop_match->push_back(false);
		     }

		 }
             */
	     }
             
	     tr.registerDerivedVec("gentop_match", gentop_match);
	     tr.registerDerivedVec("dR_top_gentop", dR_top_gentop);
	     tr.registerDerivedVec("dR_AK4_topsubjet_genmatched", dR_AK4_topsubjet_genmatched);
	     tr.registerDerivedVec("dR_AK4_top_genmatched", dR_AK4_top_genmatched);
	     tr.registerDerivedVec("top_N_AK4_matched_genmatched", top_N_AK4_matched_genmatched);
	     tr.registerDerivedVec("top_N_AK4_matched_notgenmatched", top_N_AK4_matched_notgenmatched);
	     tr.registerDerivedVec("top_N_AK4_notmatched_genmatched", top_N_AK4_notmatched_genmatched);
	     tr.registerDerivedVec("top_N_AK4_notmatched_notgenmatched", top_N_AK4_notmatched_notgenmatched);
	     tr.registerDerivedVec("top_N_AK4_matched_genmatched_0p6", top_N_AK4_matched_genmatched_0p6);
	     tr.registerDerivedVec("top_N_AK4_matched_notgenmatched_0p6", top_N_AK4_matched_notgenmatched_0p6);
	     tr.registerDerivedVec("top_N_AK4_notmatched_genmatched_0p6", top_N_AK4_notmatched_genmatched_0p6);
	     tr.registerDerivedVec("top_N_AK4_notmatched_notgenmatched_0p6", top_N_AK4_notmatched_notgenmatched_0p6);

	     tr.registerDerivedVec("top_N_AK4_matched_genmatchedother", top_N_AK4_matched_genmatchedother);
	     tr.registerDerivedVec("top_N_AK4_matched_notgenmatchedother", top_N_AK4_matched_notgenmatchedother);
	     tr.registerDerivedVec("top_N_AK4_matched_genmatchedother_0p6", top_N_AK4_matched_genmatchedother_0p6);
	     tr.registerDerivedVec("top_N_AK4_matched_notgenmatchedother_0p6", top_N_AK4_matched_notgenmatchedother_0p6);

	 }

     public:
	 Ak8DrMatch() {
	 }
	 ~Ak8DrMatch() {}
	 void operator()(NTupleReader& tr)
	 {
	     generateAk8DrMatch(tr);
	 }
     };

//void ak8DrMatch(NTupleReader& tr)
//{
   //int jetdRMatch(const std::vector<TLorentzVector>& ak8JetsLVec, const std::vector<TLorentzVector>& jetsLVec, const float jak8dRMax)
  // {
  /*
       float dRmin = 999.0;
       int minJMatch = -1;

       const int nJetsak8 = ak8JetsLVec.size();

       for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
       {
           float dR = ROOT::Math::VectorUtil::DeltaR(jetsLVec[iJet], ak8JetsLVec[iJet]);
           if(dR < dRmin)
           {
               dRmin = dR;
               minJMatch = iJet;
           }
       }
       if(dRmin < jak8dRMax) return minJMatch;
       else                return -1;
   //}
   */
//}
    //void printInterestingEvents(NTupleReader& tr)
    //{
    //    const unsigned int& run   = tr.getVar<unsigned int>("run");
    //    const unsigned int& event = tr.getVar<unsigned int>("event");
    //
    //    const float& met                            = tr.getVar<float>("met");
    //    const float& metphi                         = tr.getVar<float>("metphi");
    //
    //    const int& nMuons_CUT        = tr.getVar<int>("nMuons_CUT");
    //    const int& nElectrons_CUT    = tr.getVar<int>("nElectrons_CUT");
    //    const int& cntNJetsPt50Eta24 = tr.getVar<int>("cntNJetsPt50Eta24");
    //    const int& cntNJetsPt30Eta24 = tr.getVar<int>("cntNJetsPt30Eta24");
    //    const int& cntNJetsPt30      = tr.getVar<int>("cntNJetsPt30");
    //
    //    const float& mht    = tr.getVar<float>("mht");
    //    const float& mhtphi = tr.getVar<float>("mhtphi");
    //    const float& ht     = tr.getVar<float>("ht");
    //
    //    //if(met > 1000) std::cout << "run: " << run << "\tevent: " << event << "\tmet: " << met << "\tmetphi: " << metphi << "\tnMuons_CUT: " << nMuons_CUT << "\t nElectrons_CUT: " << nElectrons_CUT << "\tcntNJetsPt30: " << cntNJetsPt30 << "\tcntNJetsPt30Eta24: " << cntNJetsPt30Eta24 << "\tcntNJetsPt50Eta24: " << cntNJetsPt50Eta24 << "\tmht: " << mht << "\tmhtphi: " << mhtphi << "\tht: " << ht << std::endl;
    //}

}

#endif
