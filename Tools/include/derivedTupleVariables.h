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
            bool passMuTrigger     = true;
            bool passElecTrigger   = true;
            bool passMETMHTTrigger = true;
            bool passSearchTrigger = true;
            bool passHighHtTrigger = true;
            bool passPhotonTrigger = true;

            if(!tr.checkBranch("GenPart_pt"))
            {
                if( (tr.checkBranch("HLT_PFMET170_NoiseCleaned")             && tr.getVar<bool>("HLT_PFMET170_NoiseCleaned")) ||
                    (tr.checkBranch("HLT_PFMET170_JetIdCleaned")             && tr.getVar<bool>("HLT_PFMET170_JetIdCleaned")) ||
                    (tr.checkBranch("HLT_PFMET170_HBHECleaned")              && tr.getVar<bool>("HLT_PFMET170_HBHECleaned")) ||
                    (tr.checkBranch("HLT_PFMET100_PFMHT100_IDTight")         && tr.getVar<bool>("HLT_PFMET100_PFMHT100_IDTight")) ||
                    (tr.checkBranch("HLT_PFMET110_PFMHT110_IDTight")         && tr.getVar<bool>("HLT_PFMET110_PFMHT110_IDTight")) ||
                    (tr.checkBranch("HLT_PFMET120_PFMHT120_IDTight")         && tr.getVar<bool>("HLT_PFMET120_PFMHT120_IDTight")) ||
                    (tr.checkBranch("HLT_PFMET130_PFMHT130_IDTight")         && tr.getVar<bool>("HLT_PFMET130_PFMHT130_IDTight")) ||
                    (tr.checkBranch("HLT_PFMET140_PFMHT140_IDTight")         && tr.getVar<bool>("HLT_PFMET140_PFMHT140_IDTight")) ||
                    (tr.checkBranch("HLT_PFMET150_PFMHT150_IDTight")         && tr.getVar<bool>("HLT_PFMET150_PFMHT150_IDTight")) ||
                    (tr.checkBranch("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight") && tr.getVar<bool>("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight")) ||
                    (tr.checkBranch("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") && tr.getVar<bool>("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight")) ||
                    (tr.checkBranch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight") && tr.getVar<bool>("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"))
                    )
                {
                    passSearchTrigger = true;
                }

                if( (tr.checkBranch("HLT_PFHT750_4JetPt50") && tr.getVar<bool>("HLT_PFHT750_4JetPt50")) ||
                    (tr.checkBranch("HLT_PFHT800")          && tr.getVar<bool>("HLT_PFHT800")) ||
                    (tr.checkBranch("HLT_PFHT900")          && tr.getVar<bool>("HLT_PFHT900")) ||
                    (tr.checkBranch("HLT_PFJet450")         && tr.getVar<bool>("HLT_PFJet450"))
                    )
                {
                    passHighHtTrigger = true;
                }

                if( (tr.checkBranch("HLT_Photon175")                && tr.getVar<bool>("HLT_Photon175")) ||
                    (tr.checkBranch("HLT_Photon75")                 && tr.getVar<bool>("HLT_Photon75")) ||
                    (tr.checkBranch("HLT_Photon90_CaloIdL_PFHT500") && tr.getVar<bool>("HLT_Photon90_CaloIdL_PFHT500")) ||
                    (tr.checkBranch("HLT_Photon90")                 && tr.getVar<bool>("HLT_Photon90"))
                    )
                {
                    passPhotonTrigger = true;
                }

                if( (tr.checkBranch("HLT_IsoMu24")   && tr.getVar<bool>("HLT_IsoMu24")) ||
                    (tr.checkBranch("HLT_IsoTkMu24") && tr.getVar<bool>("HLT_IsoTkMu24")) ||
                    (tr.checkBranch("HLT_Mu50")      && tr.getVar<bool>("HLT_Mu50")) ||
                    (tr.checkBranch("HLT_Mu55")      && tr.getVar<bool>("HLT_Mu55"))
                    )
                {
                    passMuTrigger = true;
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
                if( tr->getVar<unsigned int>("run") >= 100000 ){ // hack to know if it's data or MC...

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
                
                std::string type;
                tr->getType("Flag_BadPFMuonFilter", type);
                bool BadPFMuonFilter = true;
                bool BadChargedCandidateFilter = true;
                if(type.compare("bool") == 0)
                {
                    BadPFMuonFilter = tr->getVar<bool>("Flag_BadPFMuonFilter");
                    BadChargedCandidateFilter = tr->getVar<bool>("Flag_BadChargedCandidateFilter");
                }
                else if(type.compare("unsigned char") == 0)
                {
                    BadPFMuonFilter = tr->getVar<unsigned char>("Flag_BadPFMuonFilter");
                    BadChargedCandidateFilter = tr->getVar<unsigned char>("Flag_BadChargedCandidateFilter");
                }
                else
                {
                    BadPFMuonFilter = tr->getVar<char>("Flag_BadPFMuonFilter");
                    BadChargedCandidateFilter = tr->getVar<char>("Flag_BadChargedCandidateFilter");
                }

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

            const std::vector<TLorentzVector>& jetsLVecNominal = tr.getVec_LVFromNano<float>("Jet");

            std::vector<TLorentzVector>& jetsLVec = tr.createDerivedVec<TLorentzVector>("jetsLVec");
            jetsLVec.reserve(jetsLVecNominal.size());
            if(tr.checkBranch("Jet_jecUncertTotal"))
            {
                const std::vector<float>& recoJetsJecUnc = tr.getVec<float>("Jet_jecUncertTotal");

                for(int ijet=0; ijet<jetsLVecNominal.size(); ++ijet)
                {
                    if(handleSys)
                    {
                        jetsLVec.push_back( jetsLVecNominal[ijet] * (1 + (JECSys * recoJetsJecUnc[ijet])));
                    }
                    else
                    {
                        jetsLVec.push_back( jetsLVecNominal[ijet] );
                    }
                }
            }
            else
            {
                for(int ijet=0; ijet<jetsLVecNominal.size(); ++ijet) jetsLVec.push_back( jetsLVecNominal[ijet] );
                handleSys = false; //If this is data, we shouldn't do anything
            }

            const std::vector<float>& recoJetsBtag      = tr.getVec<float>("Jet_btagDeepB");

            int cntCSVS = AnaFunctions::countCSVS(jetsLVec, recoJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);

            const float& met    = tr.getVar<float>("MET_pt");
            const float& metphi = tr.getVar<float>("MET_phi");

            TLorentzVector metLV;
            metLV.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

            const std::vector<TLorentzVector>& gammaLVec = tr.getVec_LVFromNano<float>("Photon");

            const std::vector<int>& photonID = tr.getVec<int>("Photon_cutBased");

            std::vector<TLorentzVector>& tightPhotons = tr.createDerivedVec<TLorentzVector>("tightPhotons");

            for(int i = 0; i < gammaLVec.size(); ++i)
            {
                if(photonID[i] == 3) // check for tight ID
                {
                    tightPhotons.push_back(gammaLVec[i]);
                }
            }

            tr.registerDerivedVar("passPhoton200", (tightPhotons.size() > 0) && (tightPhotons[0].Pt() > 200));

            const std::vector<TLorentzVector>& muonsLVec    = tr.getVec_LVFromNano<float>("Muon");
            const std::vector<float>& muonsMiniIso          = tr.getVec<float>("Muon_miniPFRelIso_all");
            const std::vector<unsigned char>& muonsMediumID         = tr.getVec<unsigned char>("Muon_mediumId");

            std::vector<float> muonsMTlep(muonsLVec.size());
            for(int i = 0; i < muonsMTlep.size(); ++i)
            {
                //muonsMTlep[i] = (muonsLVec[i] + metLV).Mt();
                muonsMTlep[i] = sqrt(2*metLV.Pt()*muonsLVec[i].Pt()*(1-cos(ROOT::Math::VectorUtil::DeltaPhi(metLV, muonsLVec[i]))));
            }

            std::vector<TLorentzVector>& cutMuVec = tr.createDerivedVec<TLorentzVector>("cutMuVec");
            std::vector<float>& cutMuMTlepVec = tr.createDerivedVec<float>("cutMuMTlepVec");
            for(int i = 0; i < muonsLVec.size(); ++i)
            {
                if(AnaFunctions::passMuon( muonsLVec[i], muonsMiniIso[i], 0.0, muonsMediumID[i], AnaConsts::muonsMiniIsoArr))
                {
                    cutMuVec.push_back(muonsLVec[i]);
                    cutMuMTlepVec.push_back(muonsMTlep[i]);
                }
            }

            const std::vector<TLorentzVector> elesLVec     = tr.getVec_LVFromNano<float>("Electron");
            const std::vector<float>& elesMiniIso          = tr.getVec<float>("Electron_miniPFRelIso_all");
            const std::vector<int>& elesCharge           = tr.getVec<int>("Electron_charge");
            
            std::vector<float>  elesMTlep(elesLVec.size());
            for(int i = 0; i < elesMTlep.size(); ++i)
            {
                //elesMTlep[i] = (elesLVec[i] + metLV).Mt();
                elesMTlep[i] = sqrt(2*metLV.Pt()*elesLVec[i].Pt()*(1-cos(ROOT::Math::VectorUtil::DeltaPhi(metLV, elesLVec[i]))));
            }

//            const std::vector<int> & elesFlagIDVec = tr.getVec<int>("Electron_cutBasedNoIso");
            //HACK
            const std::vector<int> & elesFlagIDVec = tr.getVec<int>("Electron_cutBased");

            //electron selection
            std::vector<TLorentzVector>& cutElecVec = tr.createDerivedVec<TLorentzVector>("cutElecVec");
            std::vector<float>& cutElecMTlepVec = tr.createDerivedVec<float>("cutElecMTlepVec");
            for(int i = 0; i < elesLVec.size(); ++i)
            {
                if(AnaFunctions::passElectron(elesLVec[i], elesMiniIso[i], -1, elesFlagIDVec[i] >= 3, AnaConsts::elesMiniIsoArr))
                {
                    cutElecVec.push_back(elesLVec[i]);
                    cutElecMTlepVec.push_back(elesMTlep[i]);
                }
            }

            //// Calculate number of leptons
            int nMuons = 0, nMuons_20GeV = 0, nMuons_30GeV = 0, nMuons_40GeV = 0, nMuons_50GeV = 0;
            for(auto& muon : cutMuVec)
            {
                ++nMuons;
                if(muon.Pt() > 20) ++nMuons_20GeV;
                if(muon.Pt() > 30) ++nMuons_30GeV;
                if(muon.Pt() > 40) ++nMuons_40GeV;
                if(muon.Pt() > 50) ++nMuons_50GeV;
            }
            int nElectrons = 0, nElectrons20 = 0, nElectrons30 = 0;
            for(auto& elec : cutElecVec)
            {
                ++nElectrons;
                if(elec.Pt() > 20) ++nElectrons20;
                if(elec.Pt() > 30) ++nElectrons30;
            }


            int nIsoTrks; 
            if( tr.checkBranch("IsoTrack_pt") )
            {
                nIsoTrks = 0;//AnaFunctions::countIsoTrks(tr.getVec_LVFromNano<float>("IsoTrack"), tr.getVec<float>("loose_isoTrks_iso"), tr.getVec<float>("loose_isoTrks_mtw"), tr.getVec<int>("loose_isoTrks_pdgId"));
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
            if(cutMuVec.size() >= 2)
            {
                Mmumu = (cutMuVec[0] + cutMuVec[1]).M();
                passfloatMuon = cutMuVec[0].Pt() > 30 && cutMuVec[1].Pt() > 20 && Mmumu > 81 && Mmumu < 101;
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
            float stored_weight = 1.0;
            if(tr.checkBranch("Generator_weight"))
            {
                stored_weight = tr.getVar<float>("Generator_weight");
            }

	    // Never apply this weight for data! In the old ntuple version <=3 this is "-1", in the newer ones it is "0"
	    if(stored_weight < 0) genWeight = -1.;

            //std::cout << genWeight << std::endl;
	    tr.registerDerivedVar("genWeight_sign", genWeight);

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

        int JECSys = 0;
        std::shared_ptr<TopTagger> ttMVA;
        TopCat topMatcher_;

        std::mt19937 generator;
        std::uniform_int_distribution<int> distribution;

        void prepareTopVars(NTupleReader& tr)
        {
            //We need to to handle JEC systematics here.

            bool handleSys = true;

            const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");


            const std::vector<float>& recoJetsBtag      = tr.getVec<float>("Jet_btagDeepB");
            const std::vector<float>& qgLikelihood      = tr.getVec<float>("Jet_qgl");
            
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

            std::vector<TLorentzVector> *genTops;
            std::vector<std::vector<const TLorentzVector*>> *genTopDaughters = nullptr;
            std::vector<Constituent> constituentsMVA;

            if(tr.checkBranch("GenPart_pt"))
            {
                const std::vector<TLorentzVector>& genDecayLVec = tr.getVec_LVFromNano<float>("GenPart");
                const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("GenPart_pdgId");
                const std::vector<int>& genDecayStatFlag        = tr.getVec<int>("GenPart_statusFlags");
                const std::vector<int>& genDecayMomIdxVec       = tr.getVec<int>("GenPart_genPartIdxMother");

                //prep input object (constituent) vector
                auto genMatchingInfo = ttUtility::GetTopdauGenLVecFromNano(genDecayLVec, genDecayPdgIdVec, genDecayStatFlag, genDecayMomIdxVec);
                genTops = new std::vector<TLorentzVector>(std::move(genMatchingInfo.first));
                genTopDaughters = new std::vector<std::vector<const TLorentzVector*>>(std::move(genMatchingInfo.second));
                myConstAK4Inputs = new ttUtility::ConstAK4Inputs<float>(jetsLVec, recoJetsBtag, qgLikelihood, *genTops, *genTopDaughters);

            }
            else
            {
                //no gen info is avaliable
                genTops = new std::vector<TLorentzVector>();
                
                myConstAK4Inputs = new ttUtility::ConstAK4Inputs<float>(jetsLVec, recoJetsBtag, qgLikelihood);
                
            }
                
            myConstAK4Inputs->addSupplamentalVector("qgLikelihood",                         tr.getVec<float>("Jet_qgl"));
            myConstAK4Inputs->addSupplamentalVector("qgPtD",                                tr.getVec<float>("Jet_qgptD"));
            myConstAK4Inputs->addSupplamentalVector("qgAxis1",                              tr.getVec<float>("Jet_qgAxis1"));
            myConstAK4Inputs->addSupplamentalVector("qgAxis2",                              tr.getVec<float>("Jet_qgAxis2"));
            myConstAK4Inputs->addSupplamentalVector("recoJetschargedHadronEnergyFraction",  tr.getVec<float>("Jet_chHEF"));
            myConstAK4Inputs->addSupplamentalVector("recoJetschargedEmEnergyFraction",      tr.getVec<float>("Jet_chEmEF"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsneutralEmEnergyFraction",      tr.getVec<float>("Jet_neEmEF"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsmuonEnergyFraction",           tr.getVec<float>("Jet_muEF"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsHFHadronEnergyFraction",       tr.getVec<float>("Jet_hfHadEF"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsHFEMEnergyFraction",           tr.getVec<float>("Jet_hfEMEF"));
            myConstAK4Inputs->addSupplamentalVector("recoJetsneutralEnergyFraction",        tr.getVec<float>("Jet_neHEF"));
            myConstAK4Inputs->addSupplamentalVector("PhotonEnergyFraction",                 tr.getVec<float>("Jet_phEF"));
            myConstAK4Inputs->addSupplamentalVector("ElectronEnergyFraction",               tr.getVec<float>("Jet_elEF"));
            myConstAK4Inputs->addSupplamentalVector("ChargedHadronMultiplicity",            tr.getVec<float>("Jet_chHadMult"));
            myConstAK4Inputs->addSupplamentalVector("NeutralHadronMultiplicity",            tr.getVec<float>("Jet_neHadMult"));
            myConstAK4Inputs->addSupplamentalVector("PhotonMultiplicity",                   tr.getVec<float>("Jet_phMult"));
            myConstAK4Inputs->addSupplamentalVector("ElectronMultiplicity",                 tr.getVec<float>("Jet_elMult"));
            myConstAK4Inputs->addSupplamentalVector("MuonMultiplicity",                     tr.getVec<float>("Jet_muMult"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVb",                             tr.getVec<float>("Jet_deepCSVb"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVc",                             tr.getVec<float>("Jet_deepCSVc"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVl",                             tr.getVec<float>("Jet_deepCSVudsg"));
            myConstAK4Inputs->addSupplamentalVector("DeepCSVbb",                            tr.getVec<float>("Jet_deepCSVbb"));
            myConstAK4Inputs->addSupplamentalVector("CvsL",                                 tr.getVec<float>("Jet_CvsL"));
            myConstAK4Inputs->addSupplamentalVector("qgMult",                               *convertTofloatandRegister(tr, "Jet_qgMult"));

            constituentsMVA = ttUtility::packageConstituents(*myConstAK4Inputs);
            //run tagger
            ttMVA->runTagger(constituentsMVA);

            delete myConstAK4Inputs;
            if(genTopDaughters) delete genTopDaughters;

            const TopTaggerResults& ttrMVA = ttMVA->getResults();

            const auto& candidateTops = ttrMVA.getTopCandidates();

            const auto& tops = ttrMVA.getTops();

            //get "best" top based upon on trijet mass 
            float bestTopMass = -9999.9;
            float bestTopEta = -9999.9;
            const TopObject* bestTopMassLV = nullptr;
            bool bestTopMassGenMatch = false;
            bool bestTopMassTopTag = false;
            float bestTopMassTopTagDisc = -999.9;
            
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
                    bestTopMassTopTagDisc = top.getDiscriminator();
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
            tr.registerDerivedVar("bestTopMassTopTagDisc", bestTopMassTopTagDisc);


            tr.registerDerivedVar("randomTopCand", randomTopCand?(randomTopCand->p()):(TLorentzVector()));
            tr.registerDerivedVar("randomTopCandNConst", randomTopCand?(randomTopCand->getNConstituents()):(-1));
            tr.registerDerivedVar("randomTopCandTopTag", randomTopCandTopTag);
            tr.registerDerivedVar("randomTopCandNotGenMatch", randomTopCandNotGenMatch);
        }


    public:
        PrepareTopVars(std::string taggerCfg = "TopTagger.cfg") : ttMVA(new TopTagger()), distribution(1,65000)
	{
            ttMVA->setCfgFile(taggerCfg);
	}

        ~PrepareTopVars()
        {
        }

	void operator()(NTupleReader& tr)
	{
            prepareTopVars(tr);
	}

    };
}

#endif
