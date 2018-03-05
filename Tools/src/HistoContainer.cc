#include "../include/HistoContainer.h"

//#include "../../../SusyAnaTools/Tools/NTupleReader.h"
//#include "../../../SusyAnaTools/Tools/samples.h"
//#include "../../../SusyAnaTools/Tools/SATException.h"
//#include "../../../ZInvisible/Tools/derivedTupleVariables.h"
//#include "../../../ZInvisible/Tools/baselineDef.h"
//#include "../../../ZInvisible/Tools/BTagCorrector.h"
//#include "../../../ZInvisible/Tools/TTbarCorrector.h"
//#include "../../../ZInvisible/Tools/ISRCorrector.h"
//#include "../../../ZInvisible/Tools/PileupWeights.h"
//#include "../../../ZInvisible/Tools/customize.h"
//
//#include "TopTaggerResults.h"
//#include "Constituent.h"

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include "math.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"

HistoContainer::HistoContainer(const std::string& csName = "") : csName_(csName)
{
    if(csName_.size() > 0) csName_ += "_";
    
    hMET       = bookHisto<TH1D>("MET",100,0, 1000);
    hNJets     = bookHisto<TH1D>("nJets",21,-0.5, 20.5);
    hNBJets    = bookHisto<TH1D>("nBJets",21,-0.5, 20.5);
    hNVertices = bookHisto<TH1D>("nVertices",61,-0.5, 60.5);
    hTopMass   = bookHisto<TH1D>("TopMass", 100, 0, 300);
    hTopP      = bookHisto<TH1D>("TopP", 100 , 0, 1000);
    hTopPt     = bookHisto<TH1D>("TopPt", 100, 0, 1000);
    hDiTopMass = bookHisto<TH1D>("DiTopMass", 100, 0, 1500);
    
    topPt   = bookHisto<TH1D>("topPt",   100,  0, 1000);
    topMass = bookHisto<TH1D>("topMass", 100,  0, 500);
    topEta  = bookHisto<TH1D>("topEta",  100, -5, 5);
    topCandPt   = bookHisto<TH1D>("topCandPt",   100,  0, 1000);
    topCandMass = bookHisto<TH1D>("topCandMass", 100,  0, 500);
    topCandEta  = bookHisto<TH1D>("topCandEta",  100, -5, 5);
    
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
    bestTopCandSumPt   = bookHisto<TH1D>("bestTopCandSumPt",   100,  0, 1000);
    bestTopCandSumMass = bookHisto<TH1D>("bestTopCandSumMass", 100,  0, 500);
    bestTopCandSumEta  = bookHisto<TH1D>("bestTopCandSumEta",  100, -5, 5);
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
    
    fakerateMET = bookHisto<TH1D>("fakerateMET", 100,0, 1000);
    fakerateNj  = bookHisto<TH1D>("fakerateNj", 21,-0.5, 20.5);
    fakerateNb  = bookHisto<TH1D>("fakerateNb", 21,-0.5, 20.5);
    
    fakerateMET2 = bookHisto<TH1D>("fakerateMET2", 100,0, 1000);
    fakerateNj2  = bookHisto<TH1D>("fakerateNj2", 21,-0.5, 20.5);
    fakerateNb2  = bookHisto<TH1D>("fakerateNb2", 21,-0.5, 20.5);
    
    massTemplateTop = bookHisto<TH1D>("massTemplateTop", 100,  0, 500);
    massTemplateNotTop = bookHisto<TH1D>("massTemplateBG", 100,  0, 500);
    
    allSumPt = bookHisto<TH1D>("allSumPt", 100,  0, 1000);
    bestSumPt = bookHisto<TH1D>("bestSumPt", 100,  0, 1000);
    genSumPt = bookHisto<TH1D>("genSumPt", 100,  0, 1000);
    
    topCandMassByPt = bookHisto<TH2D>("topCandMassByPt", 100,  0, 500, 100, 0, 1000);
    bestTopCandSumMassByPt = bookHisto<TH2D>("bestTopCandSumMassByPt", 100,  0, 500, 100, 0, 1000);
    bestTopCandSumMassRecoMatchByPt = bookHisto<TH2D>("bestTopCandSumMassRecoMatchByPt", 100,  0, 500, 100, 0, 1000);
    massTemplateTopByPt = bookHisto<TH2D>("massTemplateTopByPt", 100,  0, 500, 100, 0, 1000);
    massTemplateNotTopByPt = bookHisto<TH2D>("massTemplateBGByPt", 100,  0, 500, 100, 0, 1000);
    
    massTemplateGen0MatchByPt = bookHisto<TH2D>("massTemplateGen0MatchByPt", 100,  0, 500, 100, 0, 1000);
    massTemplateGen1MatchByPt = bookHisto<TH2D>("massTemplateGen1MatchByPt", 100,  0, 500, 100, 0, 1000);
    massTemplateGen2MatchByPt = bookHisto<TH2D>("massTemplateGen2MatchByPt", 100,  0, 500, 100, 0, 1000);
    massTemplateGen3MatchByPt = bookHisto<TH2D>("massTemplateGen3MatchByPt", 100,  0, 500, 100, 0, 1000);

}

void HistoContainer::setVar(const NTupleReader& tr)
{
    met_                 = &tr.getVar<double>(           "met");
    metphi_              = &tr.getVar<double>(           "metphi");    
    ht_                  = &tr.getVar<double>(           "HTTopTag");
    vtxSize_             = &tr.getVar<int>(              "vtxSize");
    cntCSVS_             = &tr.getVar<int>(              "cntCSVSTopTag");
    ttr_                 = tr.getVar<TopTaggerResults*>("ttrMVA");    
    cutMuVec_            = &tr.getVec<TLorentzVector>(   "cutMuVec");
    cutElecVec_          = &tr.getVec<TLorentzVector>(   "cutElecVec");    
    cntNJetsPt30Eta24_   = &tr.getVar<int>(              "cntNJetsPt30Eta24TopTag");    
    vTops_               = &tr.getVec<TLorentzVector>(   "vTopsNewMVA");    
    genTops_             = &tr.getVec<TLorentzVector>(   "genTops");
    genTopsRecoMatch_    = &tr.getVec<TLorentzVector>(   "vTopsGenMatchTriNewMVA");    
    vTopsNCandNewMVA_    = &tr.getVec<int>(              "vTopsNCandNewMVA");
    vTopsMatchNewMVA_    = &tr.getVec<int>(              "vTopsMatchNewMVABool");
    bestCandLV_          = &tr.getVar<TLorentzVector>(   "bestTopMassLV");
    bestTopMass_         = &tr.getVar<double>(           "bestTopMass");
    bestTopMassTopTag_   = &tr.getVar<bool>(             "bestTopMassTopTag");
    bestTopMassGenMatch_ = &tr.getVar<bool>(             "bestTopMassGenMatch");
}

void HistoContainer::fill(const NTupleReader& tr, const double& eWeight, TRandom* trand)
{
    setVar(tr);

    hMET->Fill(*met_, eWeight);
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
    
    for(const TLorentzVector& genTop : *genTopsRecoMatch_)
    {
        genTopMatchPt->Fill(genTop.Pt(), eWeight);
        genTopMatchMass->Fill(genTop.M(), eWeight);
        genTopMatchEta->Fill(genTop.Eta(), eWeight);
    }
    
    //fakerate histograms 
    for(unsigned int i = 0; i < vTopsNCandNewMVA_->size(); ++i)
    {
        if((*vTopsNCandNewMVA_)[i] == 3 && !(*vTopsMatchNewMVA_)[i])
    	{
    	    fakerateMET->Fill(*met_, eWeight);
	    fakerateNj->Fill(*cntNJetsPt30Eta24_, eWeight);
	    fakerateNb->Fill(*cntCSVS_, eWeight);
	    break;
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
    
    if(vTops_->size() > 0)
    {
    
        for(int tidx = 0; tidx < vTops_->size(); tidx++)
    	{
	    hTopMass->Fill((*vTops_)[tidx].M(),eWeight);
	    hTopP->Fill((*vTops_)[tidx].Rho(),eWeight);
	    hTopPt->Fill((*vTops_)[tidx].Perp(),eWeight);
    	}
    
        if(vTops_->size() == 2)
    	{
	    TLorentzVector diTop = (*vTops_)[0] + (*vTops_)[1];
	    hDiTopMass->Fill(diTop.M(),eWeight);
    	}
    }
    
    for(auto& top : ttr_->getTops())
    {
        massTemplateTop->Fill(top->p().M(), eWeight);
        massTemplateTopByPt->Fill(top->p().M(), top->p().Pt(), eWeight);
    
        topPt->Fill(top->p().Pt(), eWeight);
        topMass->Fill(top->p().M(), eWeight);
        topEta->Fill(top->p().Eta(), eWeight);    
    }
    
    for(auto& top : ttr_->getTops())
    {
        if(top->getNConstituents() == 3)
    	{
    	    fakerateMET2->Fill(*met_, eWeight);
	    fakerateNj2->Fill(*cntNJetsPt30Eta24_, eWeight);
	    fakerateNb2->Fill(*cntCSVS_, eWeight);
	    break;
    	}
    }
    
    //Find best candiate
                      
    //Find b jets
    std::vector<const Constituent*> bjets;
    for(const auto& constituent : ttr_->getConstituents())
    {
        if(constituent.getBTagDisc() > 0.8484)
    	{
    	    bjets.push_back(&constituent);
    	}
    }
    
    //Find lepton (here it is assumed there is exactly 1 lepton)
    TLorentzVector lepton;
    for(const auto& lep : *cutMuVec_)
    {
        if(lep.Pt() > 20)
    	{
    	    lepton = lep;
	    break;
    	}
    }
    for(const auto& lep : *cutElecVec_)
    {
        if(lep.Pt() > 20)
    	{
    	    lepton = lep;
	    break;
    	}
    }
    
    //met TLorentz vector
    TLorentzVector MET;
    MET.SetPtEtaPhiM(*met_, 0.0, *metphi_, 0);
    
    double bestSumPtVal = 99999.999;
    const TopObject* bestCand = nullptr;
    for(auto& topCand : ttr_->getTopCandidates())
    {
        switch(topCand.getGenTopMatches().size())
    	{
    	case 0:
    	    massTemplateGen0MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
	    break;
    	case 1:
    	    massTemplateGen1MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
	    break;
    	case 2:
    	    massTemplateGen2MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
	    break;
    	case 3:
    	    massTemplateGen3MatchByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
	    break;
    	}
    
        const auto& constituents = topCand.getConstituents();
        for(const auto& bjet : bjets)
    	{
    	    //Check that the b-jet is not inside the top candidate 
    	    if(std::find(constituents.begin(), constituents.end(), bjet) == constituents.end())
    	    {
	        double sumPt = (bjet->p() + MET + topCand.p() + lepton).Pt();
		allSumPt->Fill(sumPt, eWeight);
		if(topCand.getBestGenTopMatch() != nullptr) genSumPt->Fill(sumPt, eWeight);
		if(sumPt < bestSumPtVal)
    		{
    		    bestSumPtVal = sumPt;
		    bestCand = &topCand;
    	        }
    	    }
    	}
    
        topCandPt->Fill(topCand.p().Pt(), eWeight);
        topCandMass->Fill(topCand.p().M(), eWeight);
        topCandEta->Fill(topCand.p().Eta(), eWeight);
    
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
    
    if(bestCand)
    {
        bestTopCandSumPt->Fill(bestCand->p().Pt(), eWeight);
        bestTopCandSumMass->Fill(bestCand->p().M(), eWeight);
        bestTopCandSumEta->Fill(bestCand->p().Eta(), eWeight);
        bestTopCandSumMassByPt->Fill(bestCand->p().M(), bestCand->p().Pt(), eWeight);
        bestSumPt->Fill(bestSumPtVal, eWeight);
    
        for(const auto& topPtr : ttr_->getTops()) 
    	{
    	    if(topPtr == bestCand) 
    	    {
    	        bestTopCandSumMassRecoMatchByPt->Fill(bestCand->p().M(), bestCand->p().Pt(), eWeight);
		break;
    	    }
    	}
    }
    
    if(ttr_->getTopCandidates().size() > 0)
    {
        int nCand = trand->Integer(ttr_->getTopCandidates().size());
    
        const TopObject& topCand = ttr_->getTopCandidates()[nCand];
    
        randomTopCandPt->Fill(topCand.p().Pt(), eWeight);
        randomTopCandMass->Fill(topCand.p().M(), eWeight);;
        randomTopCandEta->Fill(topCand.p().Eta(), eWeight);;
              
        for(const auto& topPtr : ttr_->getTops()) 
    	{
    	    if(topPtr == &topCand) 
    	    {
    	        randomTopPt->Fill(topCand.p().Pt(), eWeight);
		randomTopMass->Fill(topCand.p().Pt(), eWeight);
		randomTopEta->Fill(topCand.p().Pt(), eWeight);
		break;
    	    }
    	}
    }        
}

void HistoContainer::save(TFile *f)
{
    f->cd();  
    for(TH1* hist : histos_) hist->Write();
}
