#ifndef HISTOCONTAINER_H
#define HISTOCONTAINER_H

#include "TopTagger/TopTagger/include/TopTaggerResults.h"

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
    TRandom* trand_;

    const double* met_;
    const double* metphi_;
    const double* ht_;
    const double* highestDisc_;
    const int* vtxSize_;
    const int* cntCSVS_;
    const TopTaggerResults* ttr_;
    const int* cntNJetsPt30Eta24_;
    const std::vector<TLorentzVector>* jets_;
    const TLorentzVector* lepton_;
    const TLorentzVector* bestCandLV_;
    const double* bestTopMass_;
    const bool* bestTopMassTopTag_;
    const bool* bestTopMassGenMatch_;
    const std::vector<TLorentzVector>* cutMuVec_;
    const std::vector<TLorentzVector>* cutElecVec_;
    const std::vector<TLorentzVector>* tightPhotonsVec_;
    const std::vector<TLorentzVector>* genTops_;

    template<typename H, typename... Args>
    H* bookHisto(const std::string& name, Args... args)
    {
        H* hptr = new H(name.c_str(), name.c_str(), args...);
        hptr->Sumw2();
        histos_.push_back(static_cast<TH1*>(hptr));
        return hptr;
    }

public:
    TH1 *jet1pT;
    TH1 *jet2pT;
    TH1 *jet3pT;
    TH1 *hMET;
    TH1 *hHT;
    TH1 *hHighestDisc;
    TH1 *hNJets;
    TH1 *hNBJets;
    TH1 *hNTops;
    TH1 *hNVertices;
    TH1 *hPhoton;
    TH1 *hMETTagged;
    TH1 *hHTTagged;
    TH1 *hNJetsTagged;
    TH1 *hNBJetsTagged;
    TH1 *hNVerticesTagged;
    TH1 *hMETTaggedGen;
    TH1 *hHTTaggedGen;
    TH1 *hHighestDiscTaggedGen;
    TH1 *hNJetsTaggedGen;
    TH1 *hNBJetsTaggedGen;
    TH1 *hNVerticesTaggedGen;
    TH1 *hMETTaggedNotGen;
    TH1 *hHTTaggedNotGen;
    TH1 *hHighestDiscTaggedNotGen;
    TH1 *hNJetsTaggedNotGen;
    TH1 *hNBJetsTaggedNotGen;
    TH1 *hNVerticesTaggedNotGen;
    TH1 *hMETTagged2;
    TH1 *hHTTagged2;
    TH1 *hHighestDiscTagged2;
    TH1 *hNJetsTagged2;
    TH1 *hNBJetsTagged2;
    TH1 *hNVerticesTagged2;
    TH1 *hPhotonTagged;
    TH1 *topPt, *topP, *topMass, *topEta, *topDisc, *topDiscGenMatch, *topDiscNotGenMatch;
    TH1 *topCandPt, *topCandMass, *topCandEta, *topCandDisc;
    TH1 *topPtGenMatch, *topPGenMatch, *topMassGenMatch, *topEtaGenMatch;
    TH1 *topCandPtGenMatch, *topCandMassGenMatch, *topCandEtaGenMatch, *topCandDiscGenMatch, *topCandDiscNotGenMatch;
    TH1 *topCandMaxDisc;
    TH1 *genTopPt, *genTopP, *genTopMass, *genTopEta;
    TH1 *genTopMatchPt, *genTopMatchMass, *genTopMatchEta;
    TH1 *genTopEvtnJet, *genTopEvtMET, *genTopEvtHT, *genTopEvtHighestDisc, *genTopEvtnVert;
    TH1 *genTopAcptEvtnJet, *genTopAcptEvtMET, *genTopAcptEvtHT, *genTopAcptEvtHighestDisc, *genTopAcptEvtnVert;
    TH1 *genTopMatchEvtnJet, *genTopMatchEvtMET, *genTopMatchEvtHT, *genTopMatchEvtHighestDisc, *genTopMatchEvtnVert;
    TH1 *bestTopCandPt, *bestTopCandMass, *bestTopCandEta;
    TH1 *bestTopCandAcptPt, *bestTopCandAcptMass, *bestTopCandAcptEta;
    TH1 *bestTopGenPt, *bestTopGenMass, *bestTopGenEta;
    TH1 *bestTopNotGenPt, *bestTopNotGenMass, *bestTopNotGenEta;
    TH1 *bestTopPt, *bestTopP, *bestTopMass, *bestTopEta;
    TH1 *randomTopCandPt,   *randomTopCandMass,   *randomTopCandEta, *randomTopCandDisc;
    TH1 *randomTopPt, *randomTopP, *randomTopMass, *randomTopEta, *randomTopDisc;
    TH2 *randomTopCandMassByPt, *randomTopMassByPt;
    TH1 *fakerateMET, *fakerateNj, *fakerateNb, *fakerateHT, *fakerateHighestDisc;
    TH1 *fakerateMET2, *fakerateNj2, *fakerateNb2, *fakerateNvert2, *fakerateHT2, *fakerateHighestDisc2;

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
    
    HistoContainer(const std::string& csName) : csName_(csName), trand_(nullptr)
    {
        trand_ = new TRandom3();
        
        jet1pT     = bookHisto<TH1D>("jet1pT",100,0,500);
        jet2pT     = bookHisto<TH1D>("jet1pT",100,0,500);
        jet3pT     = bookHisto<TH1D>("jet1pT",100,0,500);

        hMET       = bookHisto<TH1D>("MET",100,0, 1000);
        hHT        = bookHisto<TH1D>("HT",200,0, 4000);
        hHighestDisc = bookHisto<TH1D>("highestDisc",1000,-1,1);
        hNJets     = bookHisto<TH1D>("nJets",21,-0.5, 20.5);
        hNBJets    = bookHisto<TH1D>("nBJets",21,-0.5, 20.5);
        hNTops     = bookHisto<TH1D>("nTops",6,-0.5, 5.5);
        hNVertices = bookHisto<TH1D>("nVertices",61,-0.5, 60.5);
        hPhoton    = bookHisto<TH1D>("photon",100,0, 1000);

        hMETTagged       = bookHisto<TH1D>("METTagged",100,0, 1000);
        hHTTagged        = bookHisto<TH1D>("HTTagged",200,0, 4000);
        hNJetsTagged     = bookHisto<TH1D>("nJetsTagged",21,-0.5, 20.5);
        hNBJetsTagged    = bookHisto<TH1D>("nBJetsTagged",21,-0.5, 20.5);
        hNVerticesTagged = bookHisto<TH1D>("nVerticesTagged",61,-0.5, 60.5);
        hPhotonTagged    = bookHisto<TH1D>("photonTagged",100,0, 1000);

        hMETTagged2       = bookHisto<TH1D>("METTagged2",100,0, 1000);
        hHTTagged2        = bookHisto<TH1D>("HTTagged2",200,0, 4000);
        hHighestDiscTagged2 = bookHisto<TH1D>("HighestDiscTagged2",1000,-1,1);
        hNJetsTagged2     = bookHisto<TH1D>("nJetsTagged2",21,-0.5, 20.5);
        hNBJetsTagged2    = bookHisto<TH1D>("nBJetsTagged2",21,-0.5, 20.5);
        hNVerticesTagged2 = bookHisto<TH1D>("nVerticesTagged2",61,-0.5, 60.5);

        hMETTaggedNotGen       = bookHisto<TH1D>("METTaggedNotGen",100,0, 1000);
        hHTTaggedNotGen        = bookHisto<TH1D>("HTTaggedNotGen",200,0, 4000);
        hHighestDiscTaggedNotGen = bookHisto<TH1D>("HighestDiscTaggedNotGen",1000,-1,1);
        hNJetsTaggedNotGen     = bookHisto<TH1D>("nJetsTaggedNotGen",21,-0.5, 20.5);
        hNBJetsTaggedNotGen    = bookHisto<TH1D>("nBJetsTaggedNotGen",21,-0.5, 20.5);
        hNVerticesTaggedNotGen = bookHisto<TH1D>("nVerticesTaggedNotGen",61,-0.5, 60.5);

        hMETTaggedGen       = bookHisto<TH1D>("METTaggedGen",100,0, 1000);
        hHTTaggedGen        = bookHisto<TH1D>("HTTaggedGen",200,0, 4000);
        hHighestDiscTaggedGen = bookHisto<TH1D>("HighestDiscTaggedGen",1000,-1,1);
        hNJetsTaggedGen     = bookHisto<TH1D>("nJetsTaggedGen",21,-0.5, 20.5);
        hNBJetsTaggedGen    = bookHisto<TH1D>("nBJetsTaggedGen",21,-0.5, 20.5);
        hNVerticesTaggedGen = bookHisto<TH1D>("nVerticesTaggedGen",61,-0.5, 60.5);

        genTopEvtMET     = bookHisto<TH1D>("genTopEvtMET",100,0, 1000);
        genTopEvtHT      = bookHisto<TH1D>("genTopEvtHT",200,0, 4000);
        genTopEvtHighestDisc = bookHisto<TH1D>("genTopEvtHighestDisc",1000,-1,1);
        genTopEvtnJet    = bookHisto<TH1D>("genTopEvtnJet",21,-0.5, 20.5);
        genTopEvtnVert   = bookHisto<TH1D>("genTopEvtnVert",61,-0.5, 60.5);

        genTopAcptEvtMET     = bookHisto<TH1D>("genTopAcptEvtMET",100,0, 1000);
        genTopAcptEvtHT      = bookHisto<TH1D>("genTopAcptEvtHT",200,0, 4000);
        genTopAcptEvtHighestDisc = bookHisto<TH1D>("genTopAcptEvtHighestDisc",1000,-1,1);
        genTopAcptEvtnJet    = bookHisto<TH1D>("genTopAcptEvtnJet",21,-0.5, 20.5);
        genTopAcptEvtnVert   = bookHisto<TH1D>("genTopAcptEvtnVert",61,-0.5, 60.5);

        genTopMatchEvtMET     = bookHisto<TH1D>("genTopMatchEvtMET",100,0, 1000);
        genTopMatchEvtHT      = bookHisto<TH1D>("genTopMatchEvtHT",200,0, 4000);
        genTopMatchEvtHighestDisc  = bookHisto<TH1D>("genTopMatchEvtHighestDisc",1000,-1,1);
        genTopMatchEvtnJet    = bookHisto<TH1D>("genTopMatchEvtnJet",21,-0.5, 20.5);
        genTopMatchEvtnVert   = bookHisto<TH1D>("genTopMatchEvtnVert",61,-0.5, 60.5);

        topPt   = bookHisto<TH1D>("topPt",   100,  0, 1000);
        topP    = bookHisto<TH1D>("topP",   100,  0, 1000);
        topMass = bookHisto<TH1D>("topMass", 100,  0, 500);
        topEta  = bookHisto<TH1D>("topEta",  100, -5, 5);
        topDisc  = bookHisto<TH1D>("topDisc",  1000,-1,1);
        topDiscGenMatch  = bookHisto<TH1D>("topDiscGenMatch",  1000,-1,1);
        topDiscNotGenMatch  = bookHisto<TH1D>("topDiscNotGenMatch",  1000,-1,1);
        topCandPt   = bookHisto<TH1D>("topCandPt",   100,  0, 1000);
        topCandMass = bookHisto<TH1D>("topCandMass", 100,  0, 500);
        topCandEta  = bookHisto<TH1D>("topCandEta",  100, -5, 5);
        topCandDisc = bookHisto<TH1D>("topCandDisc",  1000,-1,1);

        topPtGenMatch   = bookHisto<TH1D>("topPtGenMatch",   100,  0, 1000);
        topPGenMatch    = bookHisto<TH1D>("topPGenMatch",   100,  0, 1000);
        topMassGenMatch = bookHisto<TH1D>("topMassGenMatch", 100,  0, 500);
        topEtaGenMatch  = bookHisto<TH1D>("topEtaGenMatch",  100, -5, 5);
        topCandPtGenMatch   = bookHisto<TH1D>("topCandPtGenMatch",   100,  0, 1000);
        topCandMassGenMatch = bookHisto<TH1D>("topCandMassGenMatch", 100,  0, 500);
        topCandEtaGenMatch  = bookHisto<TH1D>("topCandEtaGenMatch",  100, -5, 5);
        topCandDiscGenMatch = bookHisto<TH1D>("topCandDiscGenMatch",  1000,-1,1);
        topCandDiscNotGenMatch = bookHisto<TH1D>("topCandDiscNotGenMatch",  1000,-1,1);
        topCandMaxDisc = bookHisto<TH1D>("topCandMaxDisc",  1000,-1,1);
        
        genTopPt   = bookHisto<TH1D>("genTopPt",   100,  0, 1000);
        genTopP    = bookHisto<TH1D>("genTopP",   100,  0, 1000);
        genTopMass = bookHisto<TH1D>("genTopMass", 100,  0, 500);
        genTopEta  = bookHisto<TH1D>("genTopEta",  100, -5, 5);
        genTopMatchPt   = bookHisto<TH1D>("genTopMatchPt",   100,  0, 1000);
        genTopMatchMass = bookHisto<TH1D>("genTopMatchMass", 100,  0, 500);
        genTopMatchEta  = bookHisto<TH1D>("genTopMatchEta",  100, -5, 5);
        
        bestTopPt   = bookHisto<TH1D>("bestTopPt",   100,  0, 1000);
        bestTopP    = bookHisto<TH1D>("bestTopP",   100,  0, 1000);
        bestTopMass = bookHisto<TH1D>("bestTopMass", 100,  0, 500);
        bestTopEta  = bookHisto<TH1D>("bestTopEta",  100, -5, 5);
        bestTopCandPt   = bookHisto<TH1D>("bestTopCandPt",   100,  0, 1000);
        bestTopCandMass = bookHisto<TH1D>("bestTopCandMass", 100,  0, 500);
        bestTopCandEta  = bookHisto<TH1D>("bestTopCandEta",  100, -5, 5);
        bestTopCandAcptPt   = bookHisto<TH1D>("bestTopCandAcptPt",   100,  0, 1000);
        bestTopCandAcptMass = bookHisto<TH1D>("bestTopCandAcptMass", 100,  0, 500);
        bestTopCandAcptEta  = bookHisto<TH1D>("bestTopCandAcptEta",  100, -5, 5);
        bestTopGenPt   = bookHisto<TH1D>("bestTopGenPt",   100,  0, 1000);
        bestTopGenMass = bookHisto<TH1D>("bestTopGenMass", 100,  0, 500);
        bestTopGenEta  = bookHisto<TH1D>("bestTopGenEta",  100, -5, 5);
        bestTopNotGenPt   = bookHisto<TH1D>("bestTopNotGenPt",   100,  0, 1000);
        bestTopNotGenMass = bookHisto<TH1D>("bestTopNotGenMass", 100,  0, 500);
        bestTopNotGenEta  = bookHisto<TH1D>("bestTopNotGenEta",  100, -5, 5);
        
        randomTopPt   = bookHisto<TH1D>("randomTopPt",   100,  0, 1000);
        randomTopP    = bookHisto<TH1D>("randomTopP",   100,  0, 1000);
        randomTopMass = bookHisto<TH1D>("randomTopMass", 100,  0, 500);
        randomTopEta  = bookHisto<TH1D>("randomTopEta",  100, -5, 5);
        randomTopDisc  = bookHisto<TH1D>("randomTopDisc",  1000,-1,1);
        randomTopCandPt   = bookHisto<TH1D>("randomTopCandPt",   100,  0, 1000);
        randomTopCandMass = bookHisto<TH1D>("randomTopCandMass", 100,  0, 500);
        randomTopCandEta  = bookHisto<TH1D>("randomTopCandEta",  100, -5, 5);
        randomTopCandDisc = bookHisto<TH1D>("randomTopCandDisc",  100,  0, 1);
        randomTopMassByPt = bookHisto<TH2D>("randomTopMassByPt", 100,  0, 500, 100,0,1000);
        randomTopCandMassByPt = bookHisto<TH2D>("randomTopCandMassByPt", 100,  0, 500, 100, 0, 1000);
        
        fakerateMET = bookHisto<TH1D>("fakerateMET", 100,0, 1000);
        fakerateNj  = bookHisto<TH1D>("fakerateNj",  21,-0.5, 20.5);
        fakerateNb  = bookHisto<TH1D>("fakerateNb",  21,-0.5, 20.5);
 	fakerateHT  = bookHisto<TH1D>("fakerateHT",  200,0, 4000);
 	fakerateHighestDisc  = bookHisto<TH1D>("fakerateHighestDisc",  1000,-1,1);

        fakerateMET2 = bookHisto<TH1D>("fakerateMET2", 100,0, 1000);
        fakerateNj2  = bookHisto<TH1D>("fakerateNj2",  21,-0.5, 20.5);
        fakerateNb2  = bookHisto<TH1D>("fakerateNb2",  21,-0.5, 20.5);
        fakerateHT2  = bookHisto<TH1D>("fakerateHT2",  200,0, 4000);
        fakerateHighestDisc2  = bookHisto<TH1D>("fakerateHighestDisc2",  1000,-1,1);

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

    ~HistoContainer()
    {
        delete trand_;
    }

    void setStopVar(const TUPLECLASS& tr)
    {
        met_                 = &tr.template getVar<double>("met");
        metphi_              = &tr.template getVar<double>("metphi");    
        ht_                  = &tr.template getVar<double>("HT");
        highestDisc_          = &tr.template getVar<double>("highestDisc");
        vtxSize_             = &tr.template getVar<int>("vtxSize");
        cntCSVS_             = &tr.template getVar<int>("cntCSVS");
        ttr_                 =  tr.template getVar<TopTaggerResults*>("ttrMVA"); 
        cntNJetsPt30Eta24_   = &tr.template getVar<int>("cntNJetsPt30Eta24");
        jets_                = &tr.template getVec<TLorentzVector>("jetsLVec");
        lepton_              = &tr.template getVar<TLorentzVector>("lepton");
        bestCandLV_          = &tr.template getVar<TLorentzVector>("bestTopMassLV");
        bestTopMass_         = &tr.template getVar<double>("bestTopMass");
        bestTopMassTopTag_   = &tr.template getVar<bool>("bestTopMassTopTag");
        bestTopMassGenMatch_ = &tr.template getVar<bool>("bestTopMassGenMatch");

        //cutMuVec_            = &tr.template getVec<TLorentzVector>("cutMuVec");
        //cutElecVec_          = &tr.template getVec<TLorentzVector>("cutElecVec");    
        tightPhotonsVec_     = &tr.template getVec<TLorentzVector>("tightPhotons");  
        genTops_             = &tr.template getVec<TLorentzVector>("genTops");
    }

    void setStealthVar(const TUPLECLASS& tr)
    {
        met_                 = &tr.template getVar<double>("MET");
        metphi_              = &tr.template getVar<double>("METPhi");    
        ht_                  = &tr.template getVar<double>("HT");
        highestDisc_         = &tr.template getVar<double>("highestDisc");
        vtxSize_             = &tr.template getVar<int>("NVtx");
        cntCSVS_             = &tr.template getVar<int>("NBJets_pt30");
        ttr_                 =  tr.template getVar<TopTaggerResults*>("ttr"); 
        cntNJetsPt30Eta24_   = &tr.template getVar<int>("NJets_pt30");
        lepton_              = &tr.template getVar<TLorentzVector>("singleLepton");
        bestCandLV_          = &tr.template getVar<TLorentzVector>("bestTopMassLV");
        bestTopMass_         = &tr.template getVar<double>("bestTopMass");
        bestTopMassTopTag_   = &tr.template getVar<bool>("bestTopMassTopTag");
        bestTopMassGenMatch_ = &tr.template getVar<bool>("bestTopMassGenMatch");

        //cutMuVec_            = &tr.template getVec<TLorentzVector>("cutMuVec");
        //cutElecVec_          = &tr.template getVec<TLorentzVector>("cutElecVec");    
        tightPhotonsVec_     = &tr.template getVec<TLorentzVector>("tightPhotons");  
        genTops_             = &tr.template getVec<TLorentzVector>("hadtops");
    }

    void fill(const TUPLECLASS& tr, const double& eWeight, TRandom* trand)
    {
        if ( tr.checkBranch("met") )
        {
            setStopVar(tr);
        }
        else if( tr.checkBranch("MET") )
        {
            setStealthVar(tr);
        }

        runFill(eWeight, trand);
    }

    void fill(const TUPLECLASS& tr, const double& eWeight)
    {
        if ( tr.checkBranch("met") )
        {
            setStopVar(tr);
        }
        else if( tr.checkBranch("MET") )
        {
            setStealthVar(tr);
        }

        runFill(eWeight, trand_);
    }

    void runFill(const double& eWeight, TRandom* trand)
    {
        jet1pT->Fill((*jets_)[0].Pt(), eWeight);
        jet2pT->Fill((*jets_)[1].Pt(), eWeight);
        jet3pT->Fill((*jets_)[2].Pt(), eWeight);
    
        hMET->Fill(*met_, eWeight);
        hHT->Fill(*ht_, eWeight);
        hHighestDisc->Fill(*highestDisc_, eWeight);
        hNJets->Fill(*cntNJetsPt30Eta24_, eWeight);
        hNBJets->Fill(*cntCSVS_, eWeight);
        hNTops->Fill(ttr_->getTops().size(), eWeight);
        hNVertices->Fill(*vtxSize_,eWeight);

        if((*tightPhotonsVec_).size() > 0) hPhoton->Fill((*tightPhotonsVec_)[0].Pt(),eWeight);

        //We need to know if the identified top is gen matched, this is only useful in MC
        bool genMatched = false;
        for(const auto& top : ttr_->getTops())
        {
            const auto* genTop = top->getBestGenTopMatch();
            if(genTop)
            {
                genMatched = true;
            }
        }

        if(ttr_->getTops().size() > 0){

            hMETTagged->Fill(*met_, eWeight);
            hHTTagged->Fill(*ht_, eWeight);
            hNJetsTagged->Fill(*cntNJetsPt30Eta24_, eWeight);
            hNBJetsTagged->Fill(*cntCSVS_, eWeight);
            hNVerticesTagged->Fill(*vtxSize_,eWeight);
            if((*tightPhotonsVec_).size() > 0) hPhotonTagged->Fill((*tightPhotonsVec_)[0].Pt(),eWeight);

        }

        if(ttr_->getTops().size() > 2){

            hMETTagged2->Fill(*met_, eWeight);
            hHTTagged2->Fill(*ht_, eWeight);
            hHighestDiscTagged2->Fill(*highestDisc_, eWeight);
            hNJetsTagged2->Fill(*cntNJetsPt30Eta24_, eWeight);
            hNBJetsTagged2->Fill(*cntCSVS_, eWeight);
            hNVerticesTagged2->Fill(*vtxSize_,eWeight);

        }

        if(ttr_->getTops().size() > 1 && genMatched){

            hMETTaggedGen->Fill(*met_, eWeight);
            hHTTaggedGen->Fill(*ht_, eWeight);
            hHighestDiscTaggedGen->Fill(*highestDisc_, eWeight);
            hNJetsTaggedGen->Fill(*cntNJetsPt30Eta24_, eWeight);
            hNBJetsTaggedGen->Fill(*cntCSVS_, eWeight);
            hNVerticesTaggedGen->Fill(*vtxSize_,eWeight);

        }

        if(ttr_->getTops().size() > 1 && !genMatched){

            hMETTaggedNotGen->Fill(*met_, eWeight);
            hHTTaggedNotGen->Fill(*ht_, eWeight);
            hHighestDiscTaggedNotGen->Fill(*highestDisc_, eWeight);
            hNJetsTaggedNotGen->Fill(*cntNJetsPt30Eta24_, eWeight);
            hNBJetsTaggedNotGen->Fill(*cntCSVS_, eWeight);
            hNVerticesTaggedNotGen->Fill(*vtxSize_,eWeight);

        }

        //bools for genTop event variable plots
        bool genTopEvt      = false;
        bool genTopAcptEvt  = false;
        bool genTopMatchEvt = false;

        //plots for gen efficiency 
        for(const TLorentzVector& genTop : *genTops_)
        {
            genTopEvt = true;
            genTopPt->Fill(genTop.Pt(), eWeight);
            genTopP->Fill(genTop.P(), eWeight);
            genTopMass->Fill(genTop.M(), eWeight);
            genTopEta->Fill(genTop.Eta(), eWeight);
            if(fabs(genTop.Eta()) < 2.0) genTopAcptEvt = true;
        }

        for(const auto& top : ttr_->getTops())
        {
            const auto* genTop = top->getBestGenTopMatch();
            if(genTop)
            {
                genTopMatchEvt = true;
                genTopMatchPt->Fill(genTop->Pt(), eWeight);
                genTopMatchMass->Fill(genTop->M(), eWeight);
                genTopMatchEta->Fill(genTop->Eta(), eWeight);
            }
        }

        //Fill genTop event variable plots
        if(genTopEvt){
            genTopEvtMET->Fill(*met_, eWeight);
            genTopEvtHT->Fill(*ht_, eWeight);
            genTopEvtHighestDisc->Fill(*highestDisc_, eWeight);
            genTopEvtnJet->Fill(*cntNJetsPt30Eta24_, eWeight);
            genTopEvtnVert->Fill(*vtxSize_,eWeight);
        }
        
        if(genTopAcptEvt){
            genTopAcptEvtMET->Fill(*met_, eWeight);
            genTopAcptEvtHT->Fill(*ht_, eWeight);
            genTopAcptEvtHighestDisc->Fill(*highestDisc_, eWeight);
            genTopAcptEvtnJet->Fill(*cntNJetsPt30Eta24_, eWeight);
            genTopAcptEvtnVert->Fill(*vtxSize_,eWeight);
        }
        
        if(genTopMatchEvt){
            genTopMatchEvtMET->Fill(*met_, eWeight);
            genTopMatchEvtHT->Fill(*ht_, eWeight);
            genTopMatchEvtHighestDisc->Fill(*highestDisc_, eWeight);
            genTopMatchEvtnJet->Fill(*cntNJetsPt30Eta24_, eWeight);
            genTopMatchEvtnVert->Fill(*vtxSize_,eWeight);
        }
        
        //SF plots  
        if(*bestTopMass_ > 0.0)
        {
            bestTopCandPt->Fill(bestCandLV_->Pt(), eWeight);
            bestTopCandMass->Fill(bestCandLV_->M(), eWeight);
            bestTopCandEta->Fill(bestCandLV_->Eta(), eWeight);
        
            if(fabs(bestCandLV_->Eta()) < 2.0){
                bestTopCandAcptPt->Fill(bestCandLV_->Pt(), eWeight);
                bestTopCandAcptMass->Fill(bestCandLV_->M(), eWeight);
                bestTopCandAcptEta->Fill(bestCandLV_->Eta(), eWeight);
            }        

            if(*bestTopMassTopTag_)
            {
                bestTopPt->Fill(bestCandLV_->Pt(), eWeight);
                bestTopP->Fill(bestCandLV_->P(), eWeight);
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
            topP->Fill(top->p().P(), eWeight);
            topMass->Fill(top->p().M(), eWeight);
            topEta->Fill(top->p().Eta(), eWeight);
            topDisc->Fill(top->getDiscriminator(), eWeight);

            if(top->getBestGenTopMatch() != nullptr)
            {
                topPtGenMatch->Fill(top->p().Pt(), eWeight);
                topPGenMatch->Fill(top->p().P(), eWeight);
                topMassGenMatch->Fill(top->p().M(), eWeight);
                topEtaGenMatch->Fill(top->p().Eta(), eWeight);
                topDiscGenMatch->Fill(top->getDiscriminator(), eWeight);
            }
            else
            {
                topDiscNotGenMatch->Fill(top->getDiscriminator(), eWeight);
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
                fakerateHighestDisc->Fill(*highestDisc_, eWeight);
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
                fakerateHighestDisc2->Fill(*highestDisc_, eWeight);
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
        double discMax = 0.0;
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

            if(true)//passLepCand && nBConstituents <= 1)
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
                
                double discriminator = topCand.getDiscriminator();
                if(discriminator > discMax) discMax = discriminator;
                topCandDisc->Fill(discriminator, eWeight);

                if(topCand.getBestGenTopMatch() != nullptr)
                {
                    topCandPtGenMatch->Fill(topCand.p().Pt(), eWeight);
                    topCandMassGenMatch->Fill(topCand.p().M(), eWeight);
                    topCandEtaGenMatch->Fill(topCand.p().Eta(), eWeight);
                    topCandDiscGenMatch->Fill(discriminator, eWeight);
                }
                else
                {
                    topCandDiscNotGenMatch->Fill(discriminator, eWeight);
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
        topCandMaxDisc->Fill(discMax, eWeight);

        if(randCandIndicies.size() > 0)
        {
            int nCand = trand->Integer(randCandIndicies.size());

            const TopObject& topCand = ttr_->getTopCandidates()[randCandIndicies[nCand]];

            randomTopCandPt->Fill(topCand.p().Pt(), eWeight);
            randomTopCandMass->Fill(topCand.p().M(), eWeight);
            randomTopCandEta->Fill(topCand.p().Eta(), eWeight);
            randomTopCandDisc->Fill(topCand.getDiscriminator(), eWeight);
            randomTopCandMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);;

            for(const auto& topPtr : ttr_->getTops())
            {
                if(topPtr == &topCand)
                {
                    randomTopPt->Fill(topCand.p().Pt(), eWeight);
                    randomTopP->Fill(topCand.p().P(), eWeight);
                    randomTopMass->Fill(topCand.p().M(), eWeight);
                    randomTopEta->Fill(topCand.p().Eta(), eWeight);
                    randomTopDisc->Fill(topCand.getDiscriminator(), eWeight);
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
