#ifndef HISTOCONTAINER_H
#define HISTOCONTAINER_H

#include "TopTagger/TopTagger/interface/TopTaggerResults.h"

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

    const float* met_;
    const float* metphi_;
    const float* ht_;
    const float* highestDisc_;
    const int* vtxSize_;
    const int* cntCSVS_;
    const TopTaggerResults* ttr_;
    const int* cntNJetsPt30Eta24_;
    const std::vector<TLorentzVector>* jets_;
    const TLorentzVector* lepton_;
    const TLorentzVector* bestCandLV_;
    const float* bestTopMass_;
    const bool* bestTopMassTopTag_;
    const bool* bestTopMassGenMatch_;
    const float* bestTopMassTopTagDisc_;
    const std::vector<TLorentzVector>* cutMuVec_;
    const std::vector<TLorentzVector>* cutElecVec_;
    const std::vector<TLorentzVector>* tightPhotonsVec_;
    const std::vector<TLorentzVector>* genTops_;

    std::vector<float> workingPoints_;

    template<typename H, typename... Args>
    H* bookHisto(const std::string& name, Args... args)
    {
        H* hptr = new H(name.c_str(), name.c_str(), args...);
        hptr->Sumw2();
        histos_.push_back(static_cast<TH1*>(hptr));
        return hptr;
    }

public:
    TH1 *cutFlow_, *cutFlowNoWgt_, *passCuts_, *passCutsNoWgt_;
    TH1 *jet1pT, *jet2pT, *jet3pT;
    TH1 *hMET, *hHT, *hHighestDisc, *hNJets, *hNBJets, *hNTops, *hNVertices, *hPhoton;
    TH1 *genTopMatchPt,     *genTopMatchMass,   *genTopMatchEta;
    TH1 *topDiscGenMatch,   *topDiscNotGenMatch;
    TH1 *randomTopCandPt,   *randomTopCandMass, *randomTopCandEta, *randomTopCandDisc;
    TH2 *randomTopCandMassByPt;
    TH1 *genTopPt, *genTopP,  *genTopMass, *genTopEta;    
    TH1 *genTopEvtMET,        *genTopEvtHT,      *genTopEvtHighestDisc,      *genTopEvtnJet,      *genTopEvtnVert; 
    TH1 *genTopAcptEvtMET,    *genTopAcptEvtHT,  *genTopAcptEvtHighestDisc,  *genTopAcptEvtnJet,  *genTopAcptEvtnVert;
    TH1 *genTopMatchEvtMET,   *genTopMatchEvtHT, *genTopMatchEvtHighestDisc, *genTopMatchEvtnJet, *genTopMatchEvtnVert;
    TH1 *bestTopCandPt,       *bestTopCandMass,     *bestTopCandEta;
    TH1 *bestTopCandAcptPt,   *bestTopCandAcptMass, *bestTopCandAcptEta, *bestTopCandDisc;
    TH1 *hdPhiMin, *hdPhiMax, *hdPhiMinGenMatch, *hdPhiMaxGenMatch;
    TH1 *topCandPt,           *topCandMass,         *topCandEta,         *topCandDisc;
    TH1 *topCandPtGenMatch,   *topCandMassGenMatch, *topCandEtaGenMatch, *topCandDiscGenMatch, *topCandDiscNotGenMatch; 
    TH2 *topCandMassByPt;
    TH1 *topCandMaxDisc, *topCandMaxGenMatchDisc;
   
    std::vector<TH1*> topPt, topP, topMass, topEta, topDisc;
    std::vector<TH1*> topPtGenMatch, topPGenMatch, topMassGenMatch, topEtaGenMatch, topDiscGenMatchWP, topDiscNotGenMatchWP; 
    std::vector<TH1*> fakerateMET,     fakerateNj,  fakerateNb,  fakerateHT,  fakerateHighestDisc;
    std::vector<TH1*> fakerateMET2,    fakerateNj2, fakerateNb2, fakerateHT2, fakerateHighestDisc2, fakerateNvert2; 
    std::vector<TH1*> genTopMatchPtWP, genTopMatchMassWP, genTopMatchEtaWP;
    std::vector<TH1*> topCandDiscGenMatchWP, topCandDiscNotGenMatchWP;
    std::vector<TH1*> hMETTagged,       hHTTagged,       hNJetsTagged,             hNBJetsTagged,      hNVerticesTagged,    hPhotonTagged;
    std::vector<TH1*> hMETTagged2,      hHTTagged2,      hHighestDiscTagged2,      hNJetsTagged2,      hNBJetsTagged2,      hNVerticesTagged2;   
    std::vector<TH1*> hMETTaggedGen,    hHTTaggedGen,    hHighestDiscTaggedGen,    hNJetsTaggedGen,    hNBJetsTaggedGen,    hNVerticesTaggedGen;
    std::vector<TH1*> hMETTaggedNotGen, hHTTaggedNotGen, hHighestDiscTaggedNotGen, hNJetsTaggedNotGen, hNBJetsTaggedNotGen, hNVerticesTaggedNotGen;
    std::vector<TH1*> bestTopPt,       bestTopP,          bestTopMass, bestTopEta; 
    std::vector<TH1*> bestTopGenPt,    bestTopGenMass,    bestTopGenEta;
    std::vector<TH1*> bestTopNotGenPt, bestTopNotGenMass, bestTopNotGenEta;
    std::vector<TH1*> randomTopPt, randomTopP, randomTopMass, randomTopEta, randomTopDisc;
    std::vector<TH2*> randomTopMassByPt;

    HistoContainer(const std::string& csName, std::vector<float> workingPoints = {0.82, 0.85, 0.87, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99}) : csName_(csName), trand_(nullptr), workingPoints_(workingPoints)
    {
        trand_ = new TRandom3();
        
        jet1pT                 = bookHisto<TH1D>("jet1pT",100,0,500);
        jet2pT                 = bookHisto<TH1D>("jet2pT",100,0,500);
        jet3pT                 = bookHisto<TH1D>("jet3pT",100,0,500);

        hMET                   = bookHisto<TH1D>("MET",100,0, 1000);
        hHT                    = bookHisto<TH1D>("HT",200,0, 4000);
        hHighestDisc           = bookHisto<TH1D>("highestDisc",1000,-1,1);
        hNJets                 = bookHisto<TH1D>("nJets",21,-0.5, 20.5);
        hNBJets                = bookHisto<TH1D>("nBJets",21,-0.5, 20.5);
        hNTops                 = bookHisto<TH1D>("nTops",6,-0.5, 5.5);
        hNVertices             = bookHisto<TH1D>("nVertices",61,-0.5, 60.5);
        hPhoton                = bookHisto<TH1D>("photon",100,0, 1000);

        genTopMatchPt          = bookHisto<TH1D>("genTopMatchPt",   100,  0, 1000);
        genTopMatchMass        = bookHisto<TH1D>("genTopMatchMass", 100,  0, 500);
        genTopMatchEta         = bookHisto<TH1D>("genTopMatchEta",  100, -5, 5);

        topDiscGenMatch        = bookHisto<TH1D>("topDiscGenMatch",  1000, -1, 1); //
        topDiscNotGenMatch     = bookHisto<TH1D>("topDiscNotGenMatch",  1000, -1, 1); //

        randomTopCandPt        = bookHisto<TH1D>("randomTopCandPt",   100,  0, 1000);
        randomTopCandMass      = bookHisto<TH1D>("randomTopCandMass", 100,  0, 500);
        randomTopCandEta       = bookHisto<TH1D>("randomTopCandEta",  100, -5, 5);
        randomTopCandDisc      = bookHisto<TH1D>("randomTopCandDisc",  100,  0, 1);
        randomTopCandMassByPt  = bookHisto<TH2D>("randomTopCandMassByPt", 100,  0, 500, 100, 0, 1000);

        genTopPt               = bookHisto<TH1D>("genTopPt",   100,  0, 1000);
        genTopP                = bookHisto<TH1D>("genTopP",   100,  0, 1000);
        genTopMass             = bookHisto<TH1D>("genTopMass", 100,  0, 500);
        genTopEta              = bookHisto<TH1D>("genTopEta",  100, -5, 5);

        genTopEvtMET           = bookHisto<TH1D>("genTopEvtMET",100,0, 1000);
        genTopEvtHT            = bookHisto<TH1D>("genTopEvtHT",200,0, 4000);
        genTopEvtHighestDisc   = bookHisto<TH1D>("genTopEvtHighestDisc",1000,-1,1);
        genTopEvtnJet          = bookHisto<TH1D>("genTopEvtnJet",21,-0.5, 20.5);
        genTopEvtnVert         = bookHisto<TH1D>("genTopEvtnVert",61,-0.5, 60.5);

        genTopAcptEvtMET       = bookHisto<TH1D>("genTopAcptEvtMET",100,0, 1000);
        genTopAcptEvtHT        = bookHisto<TH1D>("genTopAcptEvtHT",200,0, 4000);
        genTopAcptEvtHighestDisc = bookHisto<TH1D>("genTopAcptEvtHighestDisc",1000,-1,1);
        genTopAcptEvtnJet      = bookHisto<TH1D>("genTopAcptEvtnJet",21,-0.5, 20.5);
        genTopAcptEvtnVert     = bookHisto<TH1D>("genTopAcptEvtnVert",61,-0.5, 60.5);

        genTopMatchEvtMET      = bookHisto<TH1D>("genTopMatchEvtMET",100,0, 1000);
        genTopMatchEvtHT       = bookHisto<TH1D>("genTopMatchEvtHT",200,0, 4000);
        genTopMatchEvtHighestDisc = bookHisto<TH1D>("genTopMatchEvtHighestDisc",1000,-1,1);
        genTopMatchEvtnJet     = bookHisto<TH1D>("genTopMatchEvtnJet",21,-0.5, 20.5);
        genTopMatchEvtnVert    = bookHisto<TH1D>("genTopMatchEvtnVert",61,-0.5, 60.5);

        bestTopCandPt          = bookHisto<TH1D>("bestTopCandPt",   100,  0, 1000);
        bestTopCandMass        = bookHisto<TH1D>("bestTopCandMass", 100,  0, 500);
        bestTopCandEta         = bookHisto<TH1D>("bestTopCandEta",  100, -5, 5);
        bestTopCandDisc        = bookHisto<TH1D>("bestTopCandDisc", 100,  0, 1);
        bestTopCandAcptPt      = bookHisto<TH1D>("bestTopCandAcptPt",   100,  0, 1000);
        bestTopCandAcptMass    = bookHisto<TH1D>("bestTopCandAcptMass", 100,  0, 500);
        bestTopCandAcptEta     = bookHisto<TH1D>("bestTopCandAcptEta",  100, -5, 5);

        hdPhiMin               = bookHisto<TH1D>("dPhiMin", 100, 0, 3.1415);
        hdPhiMax               = bookHisto<TH1D>("dPhiMax", 100, 0, 3.1415);
        hdPhiMinGenMatch       = bookHisto<TH1D>("dPhiMinGenMatch", 100, 0, 3.1415);
        hdPhiMaxGenMatch       = bookHisto<TH1D>("dPhiMaxGenMatch", 100, 0, 3.1415);

        topCandPt              = bookHisto<TH1D>("topCandPt",   100,  0, 1000);
        topCandMass            = bookHisto<TH1D>("topCandMass", 100,  0, 500);
        topCandEta             = bookHisto<TH1D>("topCandEta",  100, -5, 5);
        topCandDisc            = bookHisto<TH1D>("topCandDisc",  1000, -1, 1);

        topCandPtGenMatch      = bookHisto<TH1D>("topCandPtGenMatch",   100,  0, 1000);
        topCandMassGenMatch    = bookHisto<TH1D>("topCandMassGenMatch", 100,  0, 500);
        topCandEtaGenMatch     = bookHisto<TH1D>("topCandEtaGenMatch",  100, -5, 5);
        topCandDiscGenMatch    = bookHisto<TH1D>("topCandDiscGenMatch",  1000, -1, 1);
        topCandDiscNotGenMatch = bookHisto<TH1D>("topCandDiscNotGenMatch",  1000, -1, 1);
        topCandMaxDisc         = bookHisto<TH1D>("topCandMaxDisc",  1000, -1, 1);
        topCandMaxGenMatchDisc = bookHisto<TH1D>("topCandMaxGenMatchDisc",  1000, -1, 1);
        
        topCandMassByPt        = bookHisto<TH2D>("topCandMassByPt",     100,  0, 500, 100, 0, 1000);

        cutFlow_       = nullptr;
        cutFlowNoWgt_  = nullptr;
        passCuts_      = nullptr;
        passCutsNoWgt_ = nullptr;

        for(auto& wp : workingPoints_)
        {
            char wpStrC[32];
            sprintf(wpStrC, "%0.3f", wp);
            std::string wpStr(wpStrC);

            topPt                   .push_back( bookHisto<TH1D>("topPt_" + wpStr ,   100,  0, 1000));
            topP                    .push_back( bookHisto<TH1D>("topP_" + wpStr ,   100,  0, 1000));
            topMass                 .push_back( bookHisto<TH1D>("topMass_" + wpStr , 100,  0, 500));
            topEta                  .push_back( bookHisto<TH1D>("topEta_" + wpStr ,  100, -5, 5));
            topDisc                 .push_back( bookHisto<TH1D>("topDisc_" + wpStr ,  1000, -1, 1));

            topPtGenMatch           .push_back( bookHisto<TH1D>("topPtGenMatch_" + wpStr ,   100,  0, 1000));
            topPGenMatch            .push_back( bookHisto<TH1D>("topPGenMatch_" + wpStr ,   100,  0, 1000));
            topMassGenMatch         .push_back( bookHisto<TH1D>("topMassGenMatch_" + wpStr , 100,  0, 500));
            topEtaGenMatch          .push_back( bookHisto<TH1D>("topEtaGenMatch_" + wpStr ,  100, -5, 5));
            topDiscGenMatchWP       .push_back( bookHisto<TH1D>("topDiscGenMatchWP_" + wpStr ,  1000, -1, 1));
            topDiscNotGenMatchWP    .push_back( bookHisto<TH1D>("topDiscNotGenMatchWP_" + wpStr ,  1000, -1, 1));

            fakerateMET             .push_back( bookHisto<TH1D>("fakerateMET_" + wpStr , 100,0, 1000));
            fakerateNj              .push_back( bookHisto<TH1D>("fakerateNj_" + wpStr ,  21,-0.5, 20.5));
            fakerateNb              .push_back( bookHisto<TH1D>("fakerateNb_" + wpStr ,  21,-0.5, 20.5));
            fakerateHT              .push_back( bookHisto<TH1D>("fakerateHT_" + wpStr ,  200,0, 4000));
            fakerateHighestDisc     .push_back( bookHisto<TH1D>("fakerateHighestDisc_" + wpStr ,  1000,-1,1));

            fakerateMET2            .push_back( bookHisto<TH1D>("fakerateMET2_" + wpStr , 100,0, 1000));
            fakerateNj2             .push_back( bookHisto<TH1D>("fakerateNj2_" + wpStr ,  21,-0.5, 20.5));
            fakerateNb2             .push_back( bookHisto<TH1D>("fakerateNb2_" + wpStr ,  21,-0.5, 20.5));
            fakerateHT2             .push_back( bookHisto<TH1D>("fakerateHT2_" + wpStr ,  200,0, 4000));
            fakerateHighestDisc2    .push_back( bookHisto<TH1D>("fakerateHighestDisc2_" + wpStr ,  1000,-1,1));

            genTopMatchPtWP         .push_back( bookHisto<TH1D>("genTopMatchPtWP_" + wpStr ,   100,  0, 1000));
            genTopMatchMassWP       .push_back( bookHisto<TH1D>("genTopMatchMassWP_" + wpStr , 100,  0, 500));
            genTopMatchEtaWP        .push_back( bookHisto<TH1D>("genTopMatchEtaWP_" + wpStr ,  100, -5, 5));

            topCandDiscGenMatchWP   .push_back( bookHisto<TH1D>("topCandDiscGenMatchWP_" + wpStr , 1000, -1, 1));
            topCandDiscNotGenMatchWP.push_back( bookHisto<TH1D>("topCandDiscNotGenMatchWP_" + wpStr , 1000, -1, 1));

            hMETTagged              .push_back( bookHisto<TH1D>("METTagged_" + wpStr ,100,0, 1000));
            hHTTagged               .push_back( bookHisto<TH1D>("HTTagged_" + wpStr ,200,0, 4000));
            hNJetsTagged            .push_back( bookHisto<TH1D>("nJetsTagged_" + wpStr ,21,-0.5, 20.5));
            hNBJetsTagged           .push_back( bookHisto<TH1D>("nBJetsTagged_" + wpStr ,21,-0.5, 20.5));
            hNVerticesTagged        .push_back( bookHisto<TH1D>("nVerticesTagged_" + wpStr ,61,-0.5, 60.5));
            hPhotonTagged           .push_back( bookHisto<TH1D>("photonTagged_" + wpStr ,100,0, 1000));

            hMETTagged2             .push_back( bookHisto<TH1D>("METTagged2Top_" + wpStr ,100,0, 1000));
            hHTTagged2              .push_back( bookHisto<TH1D>("HTTagged2Top_" + wpStr ,200,0, 4000));
            hHighestDiscTagged2     .push_back( bookHisto<TH1D>("HighestDiscTagged2Top_" + wpStr ,1000,-1,1));
            hNJetsTagged2           .push_back( bookHisto<TH1D>("nJetsTagged2Top_" + wpStr ,21,-0.5, 20.5));
            hNBJetsTagged2          .push_back( bookHisto<TH1D>("nBJetsTagged2Top_" + wpStr ,21,-0.5, 20.5));
            hNVerticesTagged2       .push_back( bookHisto<TH1D>("nVerticesTagged2Top_" + wpStr ,61,-0.5, 60.5));

            hMETTaggedGen           .push_back( bookHisto<TH1D>("METTaggedGen_" + wpStr ,100,0, 1000));
            hHTTaggedGen            .push_back( bookHisto<TH1D>("HTTaggedGen_" + wpStr ,200,0, 4000));
            hHighestDiscTaggedGen   .push_back( bookHisto<TH1D>("HighestDiscTaggedGen_" + wpStr ,1000,-1,1));
            hNJetsTaggedGen         .push_back( bookHisto<TH1D>("nJetsTaggedGen_" + wpStr ,21,-0.5, 20.5));
            hNBJetsTaggedGen        .push_back( bookHisto<TH1D>("nBJetsTaggedGen_" + wpStr ,21,-0.5, 20.5));
            hNVerticesTaggedGen     .push_back( bookHisto<TH1D>("nVerticesTaggedGen_" + wpStr ,61,-0.5, 60.5));

            hMETTaggedNotGen        .push_back( bookHisto<TH1D>("METTaggedNotGen_" + wpStr ,100,0, 1000));
            hHTTaggedNotGen         .push_back( bookHisto<TH1D>("HTTaggedNotGen_" + wpStr ,200,0, 4000));
            hHighestDiscTaggedNotGen.push_back( bookHisto<TH1D>("HighestDiscTaggedNotGen_" + wpStr ,1000,-1,1));
            hNJetsTaggedNotGen      .push_back( bookHisto<TH1D>("nJetsTaggedNotGen_" + wpStr ,21,-0.5, 20.5));
            hNBJetsTaggedNotGen     .push_back( bookHisto<TH1D>("nBJetsTaggedNotGen_" + wpStr ,21,-0.5, 20.5));
            hNVerticesTaggedNotGen  .push_back( bookHisto<TH1D>("nVerticesTaggedNotGen_" + wpStr ,61,-0.5, 60.5));

            bestTopPt               .push_back( bookHisto<TH1D>("bestTopPt_" + wpStr ,   100,  0, 1000));
            bestTopP                .push_back( bookHisto<TH1D>("bestTopP_" + wpStr ,   100,  0, 1000));
            bestTopMass             .push_back( bookHisto<TH1D>("bestTopMass_" + wpStr , 100,  0, 500));
            bestTopEta              .push_back( bookHisto<TH1D>("bestTopEta_" + wpStr ,  100, -5, 5));

            bestTopGenPt            .push_back( bookHisto<TH1D>("bestTopGenPt_" + wpStr ,   100,  0, 1000));
            bestTopGenMass          .push_back( bookHisto<TH1D>("bestTopGenMass_" + wpStr , 100,  0, 500));
            bestTopGenEta           .push_back( bookHisto<TH1D>("bestTopGenEta_" + wpStr ,  100, -5, 5));
            
            bestTopNotGenPt         .push_back( bookHisto<TH1D>("bestTopNotGenPt_" + wpStr ,   100,  0, 1000));
            bestTopNotGenMass       .push_back( bookHisto<TH1D>("bestTopNotGenMass_" + wpStr , 100,  0, 500));
            bestTopNotGenEta        .push_back( bookHisto<TH1D>("bestTopNotGenEta_" + wpStr ,  100, -5, 5));

            randomTopPt             .push_back( bookHisto<TH1D>("randomTopPt_" + wpStr ,   100,  0, 1000));
            randomTopP              .push_back( bookHisto<TH1D>("randomTopP_" + wpStr ,   100,  0, 1000));
            randomTopMass           .push_back( bookHisto<TH1D>("randomTopMass_" + wpStr , 100,  0, 500));
            randomTopEta            .push_back( bookHisto<TH1D>("randomTopEta_" + wpStr ,  100, -5, 5));
            randomTopDisc           .push_back( bookHisto<TH1D>("randomTopDisc_" + wpStr ,  1000,-1,1));
            randomTopMassByPt       .push_back( bookHisto<TH2D>("randomTopMassByPt_" + wpStr , 100,  0, 500, 100,0,1000));

        }
    }

    ~HistoContainer()
    {
        delete trand_;
    }

    void setStopVar(const TUPLECLASS& tr)
    {
        met_                   = &tr.template getVar<float>("MET_pt");
        metphi_                = &tr.template getVar<float>("MET_pt");    
        ht_                    = &tr.template getVar<float>("HT");
        highestDisc_           = &tr.template getVar<float>("highestDisc");
        vtxSize_               = &tr.template getVar<int>("PV_npvs");
        cntCSVS_               = &tr.template getVar<int>("cntCSVS");
        ttr_                   =  tr.template getVar<TopTaggerResults*>("ttrMVA"); 
        cntNJetsPt30Eta24_     = &tr.template getVar<int>("cntNJetsPt30Eta24");
        jets_                  = &tr.template getVec<TLorentzVector>("jetsLVec");
        lepton_                = &tr.template getVar<TLorentzVector>("lepton");
        bestCandLV_            = &tr.template getVar<TLorentzVector>("bestTopMassLV");
        bestTopMass_           = &tr.template getVar<float>("bestTopMass");
        bestTopMassTopTag_     = &tr.template getVar<bool>("bestTopMassTopTag");
        bestTopMassGenMatch_   = &tr.template getVar<bool>("bestTopMassGenMatch");
        bestTopMassTopTagDisc_ = &tr.template getVar<float>("bestTopMassTopTagDisc");

        cutMuVec_              = &tr.template getVec<TLorentzVector>("cutMuVec");
        cutElecVec_            = &tr.template getVec<TLorentzVector>("cutElecVec");    
        tightPhotonsVec_       = &tr.template getVec<TLorentzVector>("tightPhotons");  
        genTops_               = &tr.template getVec<TLorentzVector>("genTops");
    }

    void setStealthStopVar(const TUPLECLASS& tr)
    {
        tr.registerDerivedVar("MET_float", float(tr.template getVar<double>("MET")));
        tr.registerDerivedVar("METPhi_float", float(tr.template getVar<double>("METPhi")));
        tr.registerDerivedVar("HT_float", float(tr.template getVar<double>("HT")));

        met_                   = &tr.template getVar<float>("MET_float");
        metphi_                = &tr.template getVar<float>("METPhi_float");
        ht_                    = &tr.template getVar<float>("HT_float");
        highestDisc_           = &tr.template getVar<float>("highestDisc");
        vtxSize_               = &tr.template getVar<int>("NVtx");
        cntCSVS_               = &tr.template getVar<int>("NGoodBJets_pt45");
        ttr_                   =  tr.template getVar<TopTaggerResults*>("ttr");
        cntNJetsPt30Eta24_     = &tr.template getVar<int>("NGoodJets_pt45");
        jets_                  = &tr.template getVec<TLorentzVector>("GoodJets_pt45_tlv");
        lepton_                = &tr.template getVar<TLorentzVector>("singleLepton");
        bestCandLV_            = &tr.template getVar<TLorentzVector>("bestTopMassLVCand");
        bestTopMass_           = &tr.template getVar<float>("bestTopMassCand");
        bestTopMassTopTag_     = &tr.template getVar<bool>("bestTopMassTopTagCand");
        bestTopMassGenMatch_   = &tr.template getVar<bool>("bestTopMassGenMatchCand");
        bestTopMassTopTagDisc_ = &tr.template getVar<float>("bestTopMassTopTagDisc");

        cutMuVec_              = &tr.template getVec<TLorentzVector>("Muons");
        cutElecVec_            = &tr.template getVec<TLorentzVector>("Electrons");
        tightPhotonsVec_       = &tr.template getVec<TLorentzVector>("tightPhotons");
        genTops_               = &tr.template getVec<TLorentzVector>("hadtops");
    }

    void fillWithCutFlow(const std::vector<std::pair<std::string, bool>>& cuts, const TUPLECLASS& tr, const float& eWeight, TRandom* trand)
    {
        const int EXTRACUT = 1; 

        if(cutFlow_ == nullptr)
        {
            cutFlow_       = bookHisto<TH1D>("cutFlow", cuts.size() + EXTRACUT, -0.5, cuts.size() + EXTRACUT - 0.5);
            cutFlowNoWgt_  = bookHisto<TH1D>("cutFlowNoWgt", cuts.size() + EXTRACUT, -0.5, cuts.size() + EXTRACUT - 0.5);
            passCuts_      = bookHisto<TH1D>("passCuts", cuts.size(), -0.5 + EXTRACUT, cuts.size() + EXTRACUT - 0.5);
            passCutsNoWgt_ = bookHisto<TH1D>("passCutsNoWgt", cuts.size(), -0.5 + EXTRACUT, cuts.size() + EXTRACUT - 0.5);

            cutFlow_->GetXaxis()->SetBinLabel(1, "all evt");
            cutFlowNoWgt_->GetXaxis()->SetBinLabel(1, "all evt");
            for(int iBin = 1; iBin <= cuts.size(); ++iBin)
            {
                cutFlow_->GetXaxis()->SetBinLabel(iBin + EXTRACUT, cuts[iBin - 1].first.c_str());
                cutFlowNoWgt_->GetXaxis()->SetBinLabel(iBin + EXTRACUT, cuts[iBin - 1].first.c_str());
                passCuts_->GetXaxis()->SetBinLabel(iBin, cuts[iBin - 1].first.c_str());
                passCutsNoWgt_->GetXaxis()->SetBinLabel(iBin, cuts[iBin - 1].first.c_str());
            }
        }

        cutFlow_->Fill(0.0, eWeight);
        cutFlowNoWgt_->Fill(0);

        bool passed = true;
        for(int i = 0; i < static_cast<int>(cuts.size()); ++i)
        {
            if(cuts[i].second)
            {
                passCuts_->Fill(i, eWeight);
                passCutsNoWgt_->Fill(i);
            }
            else
            {
                passed = false;
            }

            if(passed)
            {
                cutFlow_->Fill(i + EXTRACUT, eWeight);
                cutFlowNoWgt_->Fill(i + EXTRACUT);
            }
        }

        if(passed)
        {
            fill(tr, eWeight, trand);
        }
    }

    void fill(const TUPLECLASS& tr, const float& eWeight, TRandom* trand)
    {
        if(tr.checkBranch("met"))
        {
            setStopVar(tr);
        }
        else if( tr.checkBranch("MET") )
        {
            setStealthStopVar(tr);
        }
       
        runFill(eWeight, trand);
    }

    void fill(const TUPLECLASS& tr, const float& eWeight)
    {
        if(tr.checkBranch("met"))
        {
            setStopVar(tr);
        }
        else if( tr.checkBranch("MET") )
        {
            setStealthStopVar(tr);
        }

        runFill(eWeight, trand_);
    }

    // --------------------------------------------------
    // -- Fill all histograms
    // --------------------------------------------------
    void runFill(const float& eWeight, TRandom* trand)
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

        if((*tightPhotonsVec_).size() > 0) 
            hPhoton->Fill((*tightPhotonsVec_)[0].Pt(),eWeight);

        // --------------------
        // -- loop over tops 
        // --------------------
        bool genMatched     = false;
        bool genTopEvt      = false;
        bool genTopAcptEvt  = false; // gen tops in the eta acceptance of top tagger (|eta| < 2.0)
        bool genTopMatchEvt = false;

        TLorentzVector MET;
        MET.SetPtEtaPhiM(*met_, 0.0, *metphi_, 0);

        std::vector<int> randCandIndicies;
        int iCand     = 0;
        float discMax = 0.0, discMaxGenMatch = 0.0;

        for(const auto& top : ttr_->getTops())
        {
            const auto* genTop = top->getBestGenTopMatch();
            if(genTop)
            {
                genMatched = true;
                break;
            }
            
            if(genTop) 
            {
                genTopMatchEvt = true;
                genTopMatchPt->Fill(genTop->Pt(), eWeight);
                genTopMatchMass->Fill(genTop->M(), eWeight);
                genTopMatchEta->Fill(genTop->Eta(), eWeight);
            }
        
            // -- for MVA score distribution 
            float discriminator = top->getDiscriminator();
            if(top->getBestGenTopMatch() != nullptr)
            {
                topDiscGenMatch->Fill(discriminator, eWeight);
            }
            else
            {
                topDiscNotGenMatch->Fill(discriminator, eWeight);
            }
        } //  

        // -- Fill gen efficiency variables
        if(randCandIndicies.size() > 0 )
        {
            int nCand = trand->Integer(randCandIndicies.size());
            const TopObject& topCand = ttr_->getTopCandidates()[randCandIndicies[nCand]];
            randomTopCandPt->Fill(topCand.p().Pt(), eWeight);
            randomTopCandMass->Fill(topCand.p().M(), eWeight);
            randomTopCandEta->Fill(topCand.p().Eta(), eWeight);
            randomTopCandDisc->Fill(topCand.getDiscriminator(), eWeight);
            randomTopCandMassByPt->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);;
        }
       
        for(const TLorentzVector& genTop : *genTops_)
        {
            genTopEvt = true;
            genTopPt->Fill(genTop.Pt(), eWeight);
            genTopP->Fill(genTop.P(), eWeight);
            genTopMass->Fill(genTop.M(), eWeight);
            genTopEta->Fill(genTop.Eta(), eWeight);
            if(fabs(genTop.Eta()) < 2.0) genTopAcptEvt = true;
        }

        // -- Fill genTop event variables
        if(genTopEvt)
        {
            genTopEvtMET->Fill(*met_, eWeight);
            genTopEvtHT->Fill(*ht_, eWeight);
            genTopEvtHighestDisc->Fill(*highestDisc_, eWeight);
            genTopEvtnJet->Fill(*cntNJetsPt30Eta24_, eWeight);
            genTopEvtnVert->Fill(*vtxSize_,eWeight);
        }
        
        if(genTopAcptEvt)
        {
            genTopAcptEvtMET->Fill(*met_, eWeight);
            genTopAcptEvtHT->Fill(*ht_, eWeight);
            genTopAcptEvtHighestDisc->Fill(*highestDisc_, eWeight);
            genTopAcptEvtnJet->Fill(*cntNJetsPt30Eta24_, eWeight);
            genTopAcptEvtnVert->Fill(*vtxSize_,eWeight);
        }
        
        if(genTopMatchEvt)
        {
            genTopMatchEvtMET->Fill(*met_, eWeight);
            genTopMatchEvtHT->Fill(*ht_, eWeight);
            genTopMatchEvtHighestDisc->Fill(*highestDisc_, eWeight);
            genTopMatchEvtnJet->Fill(*cntNJetsPt30Eta24_, eWeight);
            genTopMatchEvtnVert->Fill(*vtxSize_,eWeight);
        }
        
        // -- SF plots 
        if(*bestTopMass_ > 0.0)
        {
            bestTopCandPt->Fill(bestCandLV_->Pt(), eWeight);
            bestTopCandMass->Fill(bestCandLV_->M(), eWeight);
            bestTopCandEta->Fill(bestCandLV_->Eta(), eWeight);
            bestTopCandDisc->Fill(*bestTopMassTopTagDisc_, eWeight);
        
            if(fabs(bestCandLV_->Eta()) < 2.0){
                bestTopCandAcptPt->Fill(bestCandLV_->Pt(), eWeight);
                bestTopCandAcptMass->Fill(bestCandLV_->M(), eWeight);
                bestTopCandAcptEta->Fill(bestCandLV_->Eta(), eWeight);
            }        
        }

        // -----------------
        // -- Find b jets   
        // -----------------      
        std::vector<const Constituent*> bjets;
        for(const Constituent& constituent : ttr_->getConstituents())
        {
            if(constituent.getBTagDisc() > 0.8484)
            {
                bjets.push_back(&constituent);
            }
        }

        // ------------------------------
        // -- loop over top candidates
        // ------------------------------
        for(auto& topCand : ttr_->getTopCandidates())
        {
            // -- delta R and Nb requirements  
            bool passLepCand   = lepton_->DeltaR(topCand.p()) > 2;
            int nBConstituents = topCand.getNBConstituents(0.8484);
            if(passLepCand && nBConstituents <= 1) 
                randCandIndicies.push_back(iCand);

            // -- dPhi cut tests 
            float dPhiMin = 999.99;
            float dPhiMax = -999.99;
            for(const auto* bjet : bjets)
            {
                float dPhi = fabs(ROOT::Math::VectorUtil::DeltaR(*lepton_ + bjet->p() + MET, topCand.p()));
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

                topCandPt->Fill(topCand.p().Pt(), eWeight);
                topCandMass->Fill(topCand.p().M(), eWeight);
                topCandEta->Fill(topCand.p().Eta(), eWeight);
                
                float discriminator = topCand.getDiscriminator();
                if(discriminator > discMax) discMax = discriminator;
                if(topCand.getBestGenTopMatch() != nullptr && discriminator > discMaxGenMatch) 
                    discMaxGenMatch = discriminator;
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
            }

            ++iCand;
        } //
        topCandMaxDisc->Fill(discMax, eWeight);
        topCandMaxGenMatchDisc->Fill(discMaxGenMatch, eWeight);

        // -----------------------------------------------
        // -- fill histograms inside working point loop
        // -----------------------------------------------
        for(unsigned int iWP = 0; iWP < workingPoints_.size(); ++iWP)
        {
            float wp = workingPoints_[iWP];
            int nTops = 0;

            // --------------------
            // -- loop over tops
            // --------------------
            for(const auto& top : ttr_->getTops())
            {
                //
                if(top->getDiscriminator() > wp) 
                    ++nTops;

                //
                if(top->getDiscriminator() > wp)
                {
                    topPt[iWP]->Fill(top->p().Pt(), eWeight);
                    topP[iWP]->Fill(top->p().P(), eWeight);
                    topMass[iWP]->Fill(top->p().M(), eWeight);
                    topEta[iWP]->Fill(top->p().Eta(), eWeight);
                    topDisc[iWP]->Fill(top->getDiscriminator(), eWeight);

                    if(top->getBestGenTopMatch() != nullptr)
                    {
                        topPtGenMatch[iWP]->Fill(top->p().Pt(), eWeight);
                        topPGenMatch[iWP]->Fill(top->p().P(), eWeight);
                        topMassGenMatch[iWP]->Fill(top->p().M(), eWeight);
                        topEtaGenMatch[iWP]->Fill(top->p().Eta(), eWeight);
                        topDiscGenMatchWP[iWP]->Fill(top->getDiscriminator(), eWeight);
                    }
                    else
                    {
                        topDiscNotGenMatchWP[iWP]->Fill(top->getDiscriminator(), eWeight);
                    }
                }

                // -- event-wise fakeaRate variables 
                if(top->getDiscriminator() > wp && top->getNConstituents() == 3 && top->getBestGenTopMatch() == nullptr)
                {
                    fakerateMET[iWP]->Fill(*met_, eWeight);
                    fakerateNj[iWP]->Fill(*cntNJetsPt30Eta24_, eWeight);
                    fakerateNb[iWP]->Fill(*cntCSVS_, eWeight);
                    fakerateHT[iWP]->Fill(*ht_, eWeight);
                    fakerateHighestDisc[iWP]->Fill(*highestDisc_, eWeight);
                    break;
                }
                
                if(top->getDiscriminator() > wp && top->getNConstituents() == 3)
                {
                    fakerateMET2[iWP]->Fill(*met_, eWeight);
                    fakerateNj2[iWP]->Fill(*cntNJetsPt30Eta24_, eWeight);
                    fakerateNb2[iWP]->Fill(*cntCSVS_, eWeight);
                    fakerateHT2[iWP]->Fill(*ht_, eWeight);
                    fakerateHighestDisc2[iWP]->Fill(*highestDisc_, eWeight);
                    break;
                }                                
                 
                // -- efficiency variables
                const auto* genTop = top->getBestGenTopMatch();
                if(genTop)
                {
                    genTopMatchEvt = true;
                    genTopMatchPtWP[iWP]->Fill(genTop->Pt(), eWeight);
                    genTopMatchMassWP[iWP]->Fill(genTop->M(), eWeight);
                    genTopMatchEtaWP[iWP]->Fill(genTop->Eta(), eWeight);
                }                
            } // 
            
            // -- for MVA score distribution
            for (auto& topCand : ttr_->getTopCandidates()) 
            {
                if(topCand.getDiscriminator() > wp)
                {
                    if(topCand.getBestGenTopMatch() != nullptr)
                    {
                        topCandDiscGenMatchWP[iWP]->Fill(topCand.getDiscriminator(), eWeight);
                    }
                    else
                    {
                        topCandDiscNotGenMatchWP[iWP]->Fill(topCand.getDiscriminator(), eWeight);
                    }
                }
            }
            
            if(nTops >= 1)
            {
                hMETTagged[iWP]->Fill(*met_, eWeight);
                hHTTagged[iWP]->Fill(*ht_, eWeight);
                hNJetsTagged[iWP]->Fill(*cntNJetsPt30Eta24_, eWeight);
                hNBJetsTagged[iWP]->Fill(*cntCSVS_, eWeight);
                hNVerticesTagged[iWP]->Fill(*vtxSize_,eWeight);
                if((*tightPhotonsVec_).size() > 0) hPhotonTagged[iWP]->Fill((*tightPhotonsVec_)[0].Pt(),eWeight);
            }

            if(nTops >= 2)
            {
                hMETTagged2[iWP]->Fill(*met_, eWeight);
                hHTTagged2[iWP]->Fill(*ht_, eWeight);
                hHighestDiscTagged2[iWP]->Fill(*highestDisc_, eWeight);
                hNJetsTagged2[iWP]->Fill(*cntNJetsPt30Eta24_, eWeight);
                hNBJetsTagged2[iWP]->Fill(*cntCSVS_, eWeight);
                hNVerticesTagged2[iWP]->Fill(*vtxSize_,eWeight);
            }

            if(nTops >= 2 && genMatched)
            {
                hMETTaggedGen[iWP]->Fill(*met_, eWeight);
                hHTTaggedGen[iWP]->Fill(*ht_, eWeight);
                hHighestDiscTaggedGen[iWP]->Fill(*highestDisc_, eWeight);
                hNJetsTaggedGen[iWP]->Fill(*cntNJetsPt30Eta24_, eWeight);
                hNBJetsTaggedGen[iWP]->Fill(*cntCSVS_, eWeight);
                hNVerticesTaggedGen[iWP]->Fill(*vtxSize_,eWeight);
            }

            if(nTops >= 2 && !genMatched)
            {
                hMETTaggedNotGen[iWP]->Fill(*met_, eWeight);
                hHTTaggedNotGen[iWP]->Fill(*ht_, eWeight);
                hHighestDiscTaggedNotGen[iWP]->Fill(*highestDisc_, eWeight);
                hNJetsTaggedNotGen[iWP]->Fill(*cntNJetsPt30Eta24_, eWeight);
                hNBJetsTaggedNotGen[iWP]->Fill(*cntCSVS_, eWeight);
                hNVerticesTaggedNotGen[iWP]->Fill(*vtxSize_,eWeight);
            }

            // -- SF plots  
            if(*bestTopMass_ > 0.0)
            {
                if(*bestTopMassTopTagDisc_ > wp)
                {
                    if(*bestTopMassTopTag_)
                    {
                        bestTopPt[iWP]->Fill(bestCandLV_->Pt(), eWeight);
                        bestTopP[iWP]->Fill(bestCandLV_->P(), eWeight);
                        bestTopMass[iWP]->Fill(bestCandLV_->M(), eWeight);
                        bestTopEta[iWP]->Fill(bestCandLV_->Eta(), eWeight);
                    }

                    if(*bestTopMassGenMatch_)
                    {
                        bestTopGenPt[iWP]->Fill(bestCandLV_->Pt(), eWeight);
                        bestTopGenMass[iWP]->Fill(bestCandLV_->M(), eWeight);
                        bestTopGenEta[iWP]->Fill(bestCandLV_->Eta(), eWeight);
                    }
                    else
                    {
                        bestTopNotGenPt[iWP]->Fill(bestCandLV_->Pt(), eWeight);
                        bestTopNotGenMass[iWP]->Fill(bestCandLV_->M(), eWeight);
                        bestTopNotGenEta[iWP]->Fill(bestCandLV_->Eta(), eWeight);
                    }
                }
            }

            // -------------------------------------
            // -- object-wise fakeaRate variables 
            // ------------------------------------- 
            if(randCandIndicies.size() > 0)
            {
                int nCand                = trand->Integer(randCandIndicies.size());
                const TopObject& topCand = ttr_->getTopCandidates()[randCandIndicies[nCand]];

                for(const auto& topPtr : ttr_->getTops())
                {
                    if(topPtr->getDiscriminator() > wp && topPtr == &topCand)
                    {
                        randomTopPt[iWP]->Fill(topCand.p().Pt(), eWeight);
                        randomTopP[iWP]->Fill(topCand.p().P(), eWeight);
                        randomTopMass[iWP]->Fill(topCand.p().M(), eWeight);
                        randomTopEta[iWP]->Fill(topCand.p().Eta(), eWeight);
                        randomTopDisc[iWP]->Fill(topCand.getDiscriminator(), eWeight);
                        randomTopMassByPt[iWP]->Fill(topCand.p().M(), topCand.p().Pt(), eWeight);
                        break;
                    }
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
