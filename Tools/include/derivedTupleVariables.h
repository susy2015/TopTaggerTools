#define TUPLE_NEW

#ifndef DERIVEDTUPLEVARIABLES_H
#define DERIVEDTUPLEVARIABLES_H

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "TopTagger/Tools/cpp/TaggerUtility.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"

#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TTModule.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"

#include "TopTagger/TopTagger/interface/TopObject.h"
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
            bool passMuTrigger = false;
            bool passElecTrigger = false;
            bool passMETMHTTrigger = false;
            bool passSearchTrigger = false, passHighHtTrigger = false, passPhotonTrigger = false;

            if( tr.getVar<bool>("HLT_PFMET170_NoiseCleaned") ||
                tr.getVar<bool>("HLT_PFMET170_JetIdCleaned") ||
                tr.getVar<bool>("HLT_PFMET170_HBHECleaned") ||
                tr.getVar<bool>("HLT_PFMET100_PFMHT100_IDTight") ||
                tr.getVar<bool>("HLT_PFMET110_PFMHT110_IDTight") ||
                tr.getVar<bool>("HLT_PFMET120_PFMHT120_IDTight") ||
                tr.getVar<bool>("HLT_PFMET130_PFMHT130_IDTight") ||
                tr.getVar<bool>("HLT_PFMET140_PFMHT140_IDTight") ||
                tr.getVar<bool>("HLT_PFMET150_PFMHT150_IDTight") ||
                tr.getVar<bool>("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight") ||
                tr.getVar<bool>("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") ||
                tr.getVar<bool>("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight")
                )
            {
                passSearchTrigger = true;
            }

            if( tr.getVar<bool>("HLT_PFHT750_4JetPt50") ||
                tr.getVar<bool>("HLT_PFHT800") ||
                tr.getVar<bool>("HLT_PFHT900") ||
                tr.getVar<bool>("HLT_PFJet450")
                )
            {
                passHighHtTrigger = true;
            }

            if( tr.getVar<bool>("HLT_Photon175") ||
                tr.getVar<bool>("HLT_Photon75") ||
                tr.getVar<bool>("HLT_Photon90_CaloIdL_PFHT500") ||
                tr.getVar<bool>("HLT_Photon90")
                )
            {
                passPhotonTrigger = true;
            }

            if( tr.getVar<bool>("HLT_IsoMu24") ||
                tr.getVar<bool>("HLT_IsoTkMu24") ||
                tr.getVar<bool>("HLT_Mu50") ||
                tr.getVar<bool>("HLT_Mu55")
                )
            {
                passMuTrigger = true;
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
            const float& met                            = tr.getVar<float>("MET_pt");
            const float& ht                             = tr.getVar<float>("HT");
            const std::vector<TLorentzVector>& cutMuVec  = tr.getVec_LVFromNano<float>("Muon");

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
                if( tr->getVar<unsigned int>("Run") >= 100000 ){ // hack to know if it's data or MC...

                    auto& goodVertexFilter = tr->getVar<bool>("Flag_goodVertices");
                    auto& globalTightHalo2016Filter = tr->getVar<bool>("Flag_globalTightHalo2016Filter");
                    auto& eeBadScFilter = tr->getVar<bool>("Flag_eeBadScFilter");

                    passDataSpec = goodVertexFilter && globalTightHalo2016Filter && eeBadScFilter;
                }

                //Check that all jets pass loose jet ID
                auto& jetID = tr->getVec<int>("Jet_jetId");
                bool passJetIDFilter = true;
                for(auto& id : jetID)
                {
                    if(!(id & 0x1)) // bit 0 is for loose ID, bit 1 for tight
                    {
                        passJetIDFilter = false;
                        break;
                    }
                }

                bool passMETratioFilter = tr->getVar<float>("CaloMET_pt")!=0 ? tr->getVar<float>("MET_pt")/tr->getVar<float>("CaloMET_pt") < 5 : true;

                auto& HBHENoiseFilter = tr->getVar<bool>("Flag_HBHENoiseFilter");
                auto& HBHENoiseIsoFilter = tr->getVar<bool>("Flag_HBHENoiseIsoFilter");
                auto& EcalDeadCellTriggerPrimitiveFilter = tr->getVar<bool>("Flag_EcalDeadCellTriggerPrimitiveFilter");
                auto& BadPFMuonFilter = tr->getVar<bool>("Flag_BadPFMuonFilter");
                auto& BadChargedCandidateFilter = tr->getVar<bool>("Flag_BadChargedCandidateFilter");

                return passDataSpec && passJetIDFilter && HBHENoiseFilter && HBHENoiseIsoFilter && EcalDeadCellTriggerPrimitiveFilter && BadPFMuonFilter && BadChargedCandidateFilter && passMETratioFilter;
            }
            catch(...)
            {
                std::cout << "PAMIC: Some filter variable is not found, ignoring event filters!!!!" << std::endl;
            }
            return true;
        }


        void prepareTopCRSelection(NTupleReader& tr)
        {
            //We need to to handle JEC systematics here.

            bool handleSys = true;

            const std::vector<TLorentzVector>& jetsLVecTemp = tr.getVec_LVFromNano<float>("Jet");
//            const std::vector<float>& recoJetsJecUnc = tr.getVec<float>("recoJetsJecUnc");

//            if(jetsLVecTemp.size() != recoJetsJecUnc.size()) handleSys = false; //If this is data, we can't do anything

            std::vector<TLorentzVector> jetsLVec;
            for(int ijet=0; ijet<jetsLVecTemp.size(); ++ijet)
            {
                if(handleSys){jetsLVec.push_back( jetsLVecTemp[ijet] * (1 + (JECSys * recoJetsJecUnc[ijet])));}
                else {jetsLVec.push_back( jetsLVecTemp[ijet] );}
            }            

            const std::vector<float>& recoJetsBtag      = tr.getVec<float>("Jet_btagCSVV2");

	    const float& stored_weight = tr.getVar<float>("genWeight");

            int cntCSVS = AnaFunctions::countCSVS(jetsLVec, recoJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);

            const float& metphi = tr.getVar<float>("MET_phi");


            const std::vector<TLorentzVector>& gammaLVec = tr.getVec_LVFromNano<float>("Photon");

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
            const std::vector<unsigned int>& elesisEB       = tr.getVec<unsigned int>("elesisEB");
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
}

#endif
