#ifndef HISTOCONTAINER_H
#define HISTOCONTAINER_H

#include "TopTaggerResults.h"

#include <vector>

#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFile.h"

template<typename TUPLECLASS>
class HistoContainer
{
private:
    std::vector<TH1*> histos_;
    std::string csName_;

    const double* met_;
    const double* metphi_;
    const double* ht_;
    const int* vtxSize_;
    const int* cntCSVS_;
    const TopTaggerResults* ttr_;
    const std::vector<TLorentzVector>* cutMuVec_;
    const std::vector<TLorentzVector>* cutElecVec_;
    const int* cntNJetsPt30Eta24_;
    const TLorentzVector* lepton_;
    const std::vector<TLorentzVector>* genTops_;
    const TLorentzVector* bestCandLV_;
    const double* bestTopMass_;
    const bool* bestTopMassTopTag_;
    const bool* bestTopMassGenMatch_;

    template<typename H, typename... Args>
    H* bookHisto(const std::string& name, Args... args)
    {
        H* hptr = new H(name.c_str(), name.c_str(), args...);
        hptr->Sumw2();
        histos_.push_back(static_cast<TH1*>(hptr));
        return hptr;
    }

public:
    TH1 *hMET;
    TH1 *hHT;
    TH1 *hNJets;
    TH1 *hNBJets;
    TH1 *hNVertices;
    TH1 *topPt, *topMass, *topEta;
    TH1 *topCandPt, *topCandMass, *topCandEta;
    TH1 *topPtGenMatch, *topMassGenMatch, *topEtaGenMatch;
    TH1 *topCandPtGenMatch, *topCandMassGenMatch, *topCandEtaGenMatch;
    TH1 *genTopPt, *genTopMass, *genTopEta;
    TH1 *genTopMatchPt, *genTopMatchMass, *genTopMatchEta;
    TH1 *bestTopCandPt, *bestTopCandMass, *bestTopCandEta;
    TH1 *bestTopGenPt, *bestTopGenMass, *bestTopGenEta;
    TH1 *bestTopNotGenPt, *bestTopNotGenMass, *bestTopNotGenEta;
    TH1 *bestTopPt, *bestTopMass, *bestTopEta;
    TH1 *randomTopCandPt,   *randomTopCandMass,   *randomTopCandEta;
    TH1 *randomTopPt, *randomTopMass, *randomTopEta;
    TH2 *randomTopCandMassByPt, *randomTopMassByPt;
    TH1 *fakerateMET, *fakerateNj, *fakerateNb, *fakerateHT;
    TH1 *fakerateMET2, *fakerateNj2, *fakerateNb2, *fakerateNvert2, *fakerateHT2;

    TH1 *massTemplateTop, *massTemplateNotTop;

    TH2 *topCandMassByPt, *massTemplateTopByPt, *massTemplateNotTopByPt;
    TH2 *massTemplateGen0MatchByPt;
    TH2 *massTemplateGen1MatchByPt;
    TH2 *massTemplateGen2MatchByPt;
    TH2 *massTemplateGen3MatchByPt;
    TH2 *massTemplateGen0MatchRecoMatchByPt;
    TH2 *massTemplateGen1MatchRecoMatchByPt;
    TH2 *massTemplateGen2MatchRecoMatchByPt;
    TH2 *massTemplateGen3MatchRecoMatchByPt;

    TH1 *hdPhiMin, *hdPhiMax;
    TH1 *hdPhiMinGenMatch, *hdPhiMaxGenMatch;
    
    HistoContainer(const std::string& csName) : csName_(csName)
    {
        hMET       = bookHisto<TH1D>("MET",100,0, 1000);
        hHT        = bookHisto<TH1D>("HT",100,0, 2000);
        hNJets     = bookHisto<TH1D>("nJets",21,-0.5, 20.5);
        hNBJets    = bookHisto<TH1D>("nBJets",21,-0.5, 20.5);
        hNVertices = bookHisto<TH1D>("nVertices",61,-0.5, 60.5);
        
        topPt   = bookHisto<TH1D>("topPt",   100,  0, 1000);
        topMass = bookHisto<TH1D>("topMass", 100,  0, 500);
        topEta  = bookHisto<TH1D>("topEta",  100, -5, 5);
        topCandPt   = bookHisto<TH1D>("topCandPt",   100,  0, 1000);
        topCandMass = bookHisto<TH1D>("topCandMass", 100,  0, 500);
        topCandEta  = bookHisto<TH1D>("topCandEta",  100, -5, 5);

        topPtGenMatch   = bookHisto<TH1D>("topPtGenMatch",   100,  0, 1000);
        topMassGenMatch = bookHisto<TH1D>("topMassGenMatch", 100,  0, 500);
        topEtaGenMatch  = bookHisto<TH1D>("topEtaGenMatch",  100, -5, 5);
        topCandPtGenMatch   = bookHisto<TH1D>("topCandPtGenMatch",   100,  0, 1000);
        topCandMassGenMatch = bookHisto<TH1D>("topCandMassGenMatch", 100,  0, 500);
        topCandEtaGenMatch  = bookHisto<TH1D>("topCandEtaGenMatch",  100, -5, 5);
        
        genTopPt   = bookHisto<TH1D>("genTopPt",   100,  0, 1000);
        genTopMass = bookHisto<TH1D>("genTopMass", 100,  0, 500);
        genTopEta  = bookHisto<TH1D>("genTopEta",  100, -5, 5);
        genTopMatchPt   = bookHisto<TH1D>("genTopMatchPt",   100,  0, 1000);
        genTopMatchMass = bookHisto<TH1D>("genTopMatchMass", 100,  0, 500);
        genTopMatchEta  = bookHisto<TH1D>("genTopMatchEta",  100, -5, 5);
        
        bestTopPt   = bookHisto<TH1D>("bestTopPt",   100,  0, 1000);
        bestTopMass = bookHisto<TH1D>("bestTopMass", 100,  0, 500);
        bestTopEta  = bookHisto<TH1D>("bestTopEta",  100, -5, 5);
        bestTopCandPt   = bookHisto<TH1D>("bestTopCandPt",   100,  0, 1000);
        bestTopCandMass = bookHisto<TH1D>("bestTopCandMass", 100,  0, 500);
        bestTopCandEta  = bookHisto<TH1D>("bestTopCandEta",  100, -5, 5);
        bestTopGenPt   = bookHisto<TH1D>("bestTopGenPt",   100,  0, 1000);
        bestTopGenMass = bookHisto<TH1D>("bestTopGenMass", 100,  0, 500);
        bestTopGenEta  = bookHisto<TH1D>("bestTopGenEta",  100, -5, 5);
        bestTopNotGenPt   = bookHisto<TH1D>("bestTopNotGenPt",   100,  0, 1000);
        bestTopNotGenMass = bookHisto<TH1D>("bestTopNotGenMass", 100,  0, 500);
        bestTopNotGenEta  = bookHisto<TH1D>("bestTopNotGenEta",  100, -5, 5);
        
        randomTopPt   = bookHisto<TH1D>("randomTopPt",   100,  0, 1000);
        randomTopMass = bookHisto<TH1D>("randomTopMass", 100,  0, 500);
        randomTopEta  = bookHisto<TH1D>("randomTopEta",  100, -5, 5);
        randomTopCandPt   = bookHisto<TH1D>("randomTopCandPt",   100,  0, 1000);
        randomTopCandMass = bookHisto<TH1D>("randomTopCandMass", 100,  0, 500);
        randomTopCandEta  = bookHisto<TH1D>("randomTopCandEta",  100, -5, 5);
        randomTopMassByPt = bookHisto<TH2D>("randomTopMassByPt", 100,  0, 500, 100, 0, 1000);
        randomTopCandMassByPt = bookHisto<TH2D>("randomTopCandMassByPt", 100,  0, 500, 100, 0, 1000);
        
        fakerateMET = bookHisto<TH1D>("fakerateMET", 100,0, 1000);
        fakerateNj  = bookHisto<TH1D>("fakerateNj",  21,-0.5, 20.5);
        fakerateNb  = bookHisto<TH1D>("fakerateNb",  21,-0.5, 20.5);
 	fakerateHT  = bookHisto<TH1D>("fakerateHT",  100,0, 2000);

        fakerateMET2 = bookHisto<TH1D>("fakerateMET2", 100,0, 1000);
        fakerateNj2  = bookHisto<TH1D>("fakerateNj2",  21,-0.5, 20.5);
        fakerateNb2  = bookHisto<TH1D>("fakerateNb2",  21,-0.5, 20.5);
        fakerateHT2  = bookHisto<TH1D>("fakerateHT2",  100,0, 2000);

        massTemplateTop = bookHisto<TH1D>("massTemplateTop", 100,  0, 500);
        massTemplateNotTop = bookHisto<TH1D>("massTemplateBG", 100,  0, 500);
        
        topCandMassByPt        = bookHisto<TH2D>("topCandMassByPt",     100,  0, 500, 100, 0, 1000);
        massTemplateTopByPt    = bookHisto<TH2D>("massTemplateTopByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateNotTopByPt = bookHisto<TH2D>("massTemplateBGByPt",  100,  0, 500, 100, 0, 1000);

        massTemplateGen0MatchByPt = bookHisto<TH2D>("massTemplateGen0MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen1MatchByPt = bookHisto<TH2D>("massTemplateGen1MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen2MatchByPt = bookHisto<TH2D>("massTemplateGen2MatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen3MatchByPt = bookHisto<TH2D>("massTemplateGen3MatchByPt", 100,  0, 500, 100, 0, 1000);
        
        massTemplateGen0MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen0MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen1MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen1MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen2MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen2MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
        massTemplateGen3MatchRecoMatchByPt = bookHisto<TH2D>("massTemplateGen3MatchRecoMatchByPt", 100,  0, 500, 100, 0, 1000);

        hdPhiMin = bookHisto<TH1D>("dPhiMin", 100, 0, 3.1415);
        hdPhiMax = bookHisto<TH1D>("dPhiMax", 100, 0, 3.1415);
        hdPhiMinGenMatch = bookHisto<TH1D>("dPhiMinGenMatch", 100, 0, 3.1415);
        hdPhiMaxGenMatch = bookHisto<TH1D>("dPhiMaxGenMatch", 100, 0, 3.1415);
    }

    void setStopVar(const TUPLECLASS& tr)
    {
        met_                 = &tr.template getVar<double>(           "met");
        metphi_              = &tr.template getVar<double>(           "metphi");    
        ht_                  = &tr.template getVar<double>(           "HT");
        vtxSize_             = &tr.template getVar<int>(              "vtxSize");
        cntCSVS_             = &tr.template getVar<int>(              "cntCSVS");
        ttr_                 =  tr.template getVar<TopTaggerResults*>("ttrMVA");    
        cntNJetsPt30Eta24_   = &tr.template getVar<int>(              "cntNJetsPt30Eta24");
        lepton_              = &tr.template getVar<TLorentzVector>(   "lepton");
        genTops_             = &tr.template getVec<TLorentzVector>(   "genTops");
        bestCandLV_          = &tr.template getVar<TLorentzVector>(   "bestTopMassLV");
        bestTopMass_         = &tr.template getVar<double>(           "bestTopMass");
        bestTopMassTopTag_   = &tr.template getVar<bool>(             "bestTopMassTopTag");
        bestTopMassGenMatch_ = &tr.template getVar<bool>(             "bestTopMassGenMatch");
    }

    void setStealthVar(const TUPLECLASS& tr)
    {
        met_                 = *tr.MET;
        metphi_              = *tr.METPhi;
        ht_                  = *tr.HT;
        vtxSize_             = *tr.NVtx;
        //cntCSVS_             = *tr.;
        //ttr_                 = *tr.;
        //cntNJetsPt30Eta24_   = *tr.;
        //lepton_              = *tr.;
        //genTops_             = *tr.;
        //bestCandLV_          = *tr.;
        //bestTopMass_         = *tr.;
        //bestTopMassTopTag_   = *tr.;
        //bestTopMassGenMatch_ = *tr.;
    }

    ///Hack
    ///Should be fix in simpleAnalyzer.C
    void fill(const TUPLECLASS& tr, const double& eWeight, TRandom* trand)
    {
        setStopVar(tr);
        runFill(eWeight, trand);
    }

    void fill(const TUPLECLASS& tr, const double& eWeight)
    {
        TRandom* trand = new TRandom3();
        setStealthVar(tr);
        runFill(eWeight, trand);
    }

    void runFill(const double& eWeight, TRandom* trand)
    {    

        hMET->Fill(*met_, eWeight);
        hHT->Fill(*ht_, eWeight);
        hNJets->Fill(*cntNJetsPt30Eta24_, eWeight);
        hNBJets->Fill(*cntCSVS_, eWeight);
        hNVertices->Fill(*vtxSize_,eWeight);

        //plots for gen efficiency 
        for(const TLorentzVector& genTop : *genTops_)
        {
            genTopPt->Fill(genTop.Pt(), eWeight);
            genTopMass->Fill(genTop.M(), eWeight);
            genTopEta->Fill(genTop.Eta(), eWeight);
        }
        
        for(const auto& top : ttr_->getTops())
        {
            const auto* genTop = top->getBestGenTopMatch();
            if(genTop)
            {
                genTopMatchPt->Fill(genTop->Pt(), eWeight);
                genTopMatchMass->Fill(genTop->M(), eWeight);
                genTopMatchEta->Fill(genTop->Eta(), eWeight);
            }
        }
        
        //SF plots  
        if(*bestTopMass_ > 0.0)
        {
            bestTopCandPt->Fill(bestCandLV_->Pt(), eWeight);
            bestTopCandMass->Fill(bestCandLV_->M(), eWeight);
            bestTopCandEta->Fill(bestCandLV_->Eta(), eWeight);

            if(*bestTopMassTopTag_)
            {
                bestTopPt->Fill(bestCandLV_->Pt(), eWeight);
                bestTopMass->Fill(bestCandLV_->M(), eWeight);
                bestTopEta->Fill(bestCandLV_->Eta(), eWeight);
            }

            if(*bestTopMassGenMatch_)
            {
                bestTopGenPt->Fill(bestCandLV_->Pt(), eWeight);
                bestTopGenMass->Fill(bestCandLV_->M(), eWeight);
                bestTopGenEta->Fill(bestCandLV_->Eta(), eWeight);
            }
            else
            {
                bestTopNotGenPt->Fill(bestCandLV_->Pt(), eWeight);
                bestTopNotGenMass->Fill(bestCandLV_->M(), eWeight);
                bestTopNotGenEta->Fill(bestCandLV_->Eta(), eWeight);
            }
        }

        for(auto& top : ttr_->getTops())
        {
            massTemplateTop->Fill(top->p().M(), eWeight);
            massTemplateTopByPt->Fill(top->p().M(), top->p().Pt(), eWeight);

            topPt->Fill(top->p().Pt(), eWeight);
            topMass->Fill(top->p().M(), eWeight);
            topEta->Fill(top->p().Eta(), eWeight);

            if(top->getBestGenTopMatch() != nullptr)
            {
                topPtGenMatch->Fill(top->p().Pt(), eWeight);
                topMassGenMatch->Fill(top->p().M(), eWeight);
                topEtaGenMatch->Fill(top->p().Eta(), eWeight);
            }
        }

        for(auto& top : ttr_->getTops())
        {
            if(top->getNConstituents() == 3 && top->getBestGenTopMatch() == nullptr)
            {
                fakerateMET->Fill(*met_, eWeight);
                fakerateNj->Fill(*cntNJetsPt30Eta24_, eWeight);
                fakerateNb->Fill(*cntCSVS_, eWeight);
                fakerateHT->Fill(*ht_, eWeight);
                break;
            }
        }

        for(auto& top : ttr_->getTops())
        {
            if(top->getNConstituents() == 3)
            {
                fakerateMET2->Fill(*met_, eWeight);
                fakerateNj2->Fill(*cntNJetsPt30Eta24_, eWeight);
                fakerateNb2->Fill(*cntCSVS_, eWeight);
                fakerateHT2->Fill(*ht_, eWeight);
                break;
            }
        }

        //Find best candiate  

        //met TLorentz vector  
        TLorentzVector MET;
        MET.SetPtEtaPhiM(*met_, 0.0, *metphi_, 0);

        //double xi = MET.Px()*lepton_->Px() + MET.Py()*lepton_->Py();
        //double mW = 80.385;
        //double a = pow(lepton_->Pt(), 2);
        //double b = -lepton_->Pz()*(mW*mW + 2*xi);
        //double c = pow(lepton_->P(), 2) * pow(*met_, 2) - pow(xi + mW*mW/2.0, 2);
        //double pnuzp = (-b + sqrt(b*b - 4*a*c))/(2*a);
        //double pnuzm = (-b - sqrt(b*b - 4*a*c))/(2*a);
        //
        //TLorentzVector neutrinop(MET.X(), MET.Y(), pnuzp, 0.0);
        //TLorentzVector neutrinom(MET.X(), MET.Y(), pnuzm, 0.0);
        //
        //std::cout << xi << "\t" << a << "\t" << b << "\t" << c << "\t" << pnuzp << "\t" << pnuzm << "\t" << (*lepton_ + neutrinop).M() << "\t" << (*lepton_ + neutrinom).M() <<  std::endl;

        //Find b jets         
        std::vector<const Constituent*> bjets;
        for(const Constituent& constituent : ttr_->getConstituents())
        {
            if(constituent.getBTagDisc() > 0.8484)
            {
                bjets.push_back(&constituent);
            }
        }


        double bestSumPtVal = 99999.999;
        const TopObject* bestCand = nullptr;
        std::vector<int> randCandIndicies;
        int iCand = 0;
        for(auto& topCand : ttr_->getTopCandidates())
        {
            //delta R and Nb requirements  
            bool passLepCand = lepton_->DeltaR(topCand.p()) > 2;
            int nBConstituents = topCand.getNBConstituents(0.8484);
            if(passLepCand && nBConstituents <= 1) randCandIndicies.push_back(iCand);

            //dPhi cut tests 
            double dPhiMin = 999.99;
            double dPhiMax = -999.99;
            
            for(const auto* bjet : bjets)
            {
                double dPhi = fabs(ROOT::Math::VectorUtil::DeltaR(*lepton_ + bjet->p() + MET, topCand.p()));
                for(const auto* constituent : topCand.getConstituents())
                {
                    if(bjet == constituent) continue;
                }
                dPhiMin = std::min(dPhiMin, dPhi);
                dPhiMax = std::max(dPhiMax, dPhi);
            }

            hdPhiMin->Fill(dPhiMin, eWeight);
            hdPhiMax->Fill(dPhiMax, eWeight);

            if(topCand.getBestGenTopMatch() != nullptr)
            {
                hdPhiMinGenMatch->Fill(dPhiMin, eWeight);
                hdPhiMaxGenMatch->Fill(dPhiMax, eWeight);
            }

            bool recoMatch = false;
            for(auto& top : ttr_->getTops())
            {
                if(top == &topCand)
                {
                    recoMatch = true;
                    break;
                }
            }

            if(passLepCand && nBConstituents <= 1)
            {
                int nGenMatch = 0;
                for(const auto& genMatch : topCand.getGenTopMatches())
                {
                    nGenMatch = std::max(nGenMatch, static_cast<int>(genMatch.second.size()));
                }
                switch(nGenMatch)
                {
                case 0:
                    massTemplateGen0MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen0MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                case 1:
                    massTemplateGen1MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen1MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                case 2:
                    massTemplateGen2MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen2MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                case 3:
                    massTemplateGen3MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    if(recoMatch) massTemplateGen3MatchRecoMatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                }

                topCandPt->Fill(topCand.p().Pt(), eWeight);
                topCandMass->Fill(topCand.p().M(), eWeight);
                topCandEta->Fill(topCand.p().Eta(), eWeight);

                if(topCand.getBestGenTopMatch() != nullptr)
                {
                    topCandPtGenMatch->Fill(topCand.p().Pt(), eWeight);
                    topCandMassGenMatch->Fill(topCand.p().M(), eWeight);
                    topCandEtaGenMatch->Fill(topCand.p().Eta(), eWeight);
                }

                topCandMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);

                if(topCand.getDiscriminator() > std::min(0.97, 0.8 + 0.0005*topCand.p().Pt()))
                {
                    //massTemplateTop->Fill(topCand.p().M(), eWeight);  
                    //massTemplateTopByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);   
                }
                else
                {
                    massTemplateNotTop->Fill(topCand.p().M(), eWeight);
                    massTemplateNotTopByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                }
            }

            ++iCand;
        }

        if(randCandIndicies.size() > 0)
        {
            int nCand = trand->Integer(randCandIndicies.size());

            const TopObject& topCand = ttr_->getTopCandidates()[randCandIndicies[nCand]];

            randomTopCandPt->Fill(topCand.p().Pt(), eWeight);
            randomTopCandMass->Fill(topCand.p().M(), eWeight);
            randomTopCandEta->Fill(topCand.p().Eta(), eWeight);
            randomTopCandMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);;

            for(const auto& topPtr : ttr_->getTops())
            {
                if(topPtr == &topCand)
                {
                    randomTopPt->Fill(topCand.p().Pt(), eWeight);
                    randomTopMass->Fill(topCand.p().M(), eWeight);
                    randomTopEta->Fill(topCand.p().Eta(), eWeight);
                    randomTopMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                    break;
                }
            }
        }
    }
    
    void save(TFile *f)
    {
        f->cd();
        TDirectory *td = f->mkdir(csName_.c_str(), csName_.c_str());
        td->cd();
        for(TH1* hist : histos_) hist->Write();
        f->cd();
    }
    
};

#endif
