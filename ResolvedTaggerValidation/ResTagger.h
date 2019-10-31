#ifndef RESTAGGER_H
#define RESTAGGER_H

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "TopTagger/Tools/cpp/TaggerUtility.h"
#include "TopTagger/Tools/cpp/PlotUtility.h"

#include "TopTagger/TopTagger/interface/TopObject.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/TTModule.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "TopTagger/CfgParser/include/CfgDocument.hh"


#include "TFile.h"
using namespace std;

 class PrepareTopVars
 {
 private:

   std::shared_ptr<TopTagger> ttMVA;
   TopCat topMatcher_;
   
   void prepareTopVars(NTupleReader& tr)
   {
     
     const std::vector<TLorentzVector>& jetsLVec = tr.getVec_LVFromNano<float>("Jet");
     const std::vector<float>& recoJetsBtag      = tr.getVec<float>("Jet_btagDeepB");
     const std::vector<float>& qgLikelihood = tr.getVec<float>("Jet_qgl");  
     
  /* //For Sys
  const std::vector<TLorentzVector>& jetsLVecTemp = tr.getVec<TLorentzVector>("jetsLVec");
  const std::vector<float>& recoJetsJecUnc = tr.getVec<float>("recoJetsJecUnc");
  std::vector<TLorentzVector> jetsLVec;
  for(int ijet=0; ijet<jetsLVecTemp.size(); ++ijet)
    {
      //Up
      jetsLVec.push_back(jetsLVecTemp[ijet] * (1 + recoJetsJecUnc[ijet]));
      //Down 
      jetsLVec.push_back(jetsLVecTemp[ijet] * (1 - recoJetsJecUnc[ijet]));
    }
  */ 
  //Helper function to turn int vectors into float vectors
  auto convertTofloatandRegister1 = [](NTupleReader& tr, const std::string& name)
    {
      const std::vector<int>& intVec = tr.getVec<int>(name);
      std::vector<float>* floatVec = new std::vector<float>(intVec.begin(), intVec.end());
      tr.registerDerivedVec(name+"ConvertedTofloat", floatVec);
      return floatVec;
    };

  //Helper function to turn double vectors into float vectors
  auto convertTofloatandRegister2 = [](NTupleReader& tr, const std::string& name)
    {
      const std::vector<double>& doubleVec = tr.getVec<double>(name);
      std::vector<float>* floatVec = new std::vector<float>(doubleVec.begin(), doubleVec.end());
      tr.registerDerivedVec(name+"ConvertedTofloat", floatVec);
      return floatVec;
    };

  //Helper function to turn int vectors into double vectors
  auto convertTodoubleandRegister1 = [](NTupleReader& tr, const std::string& name)
    {
      const std::vector<int>& intVec = tr.getVec<int>(name);
      std::vector<double>* doubleVec = new std::vector<double>(intVec.begin(), intVec.end());
      tr.registerDerivedVec(name+"ConvertedTodouble", doubleVec);
      return doubleVec;
    };

 //Helper function to turn float vectors into double vectors
  auto convertTodoubleandRegister2 = [](NTupleReader& tr, const std::string& name)
    {
      const std::vector<float>& floatVec = tr.getVec<float>(name);
      std::vector<double>* doubleVec = new std::vector<double>(floatVec.begin(), floatVec.end());
      tr.registerDerivedVec(name+"ConvertedTodouble", doubleVec);
      return doubleVec;
    };
  
  //Helper function to turn char vectors into int vectors
  auto convertTointandRegister1 = [](NTupleReader& tr, const std::string& name)
    {
      const std::vector<unsigned char>& charVec = tr.getVec<unsigned char>(name);
      std::vector<int>* intVec = new std::vector<int>(charVec.begin(), charVec.end());
      tr.registerDerivedVec(name+"ConvertedToint", intVec);
      return intVec;
    };

  //Helper function to turn int vectors into char vectors
  auto convertTocharandRegister1 = [](NTupleReader& tr, const std::string& name)
    {
      const std::vector<int>& intVec = tr.getVec<int>(name);
      std::vector<unsigned char>* charVec = new std::vector<unsigned char>(intVec.begin(), intVec.end());
      tr.registerDerivedVec(name+"ConvertedTochar", charVec);
      return charVec;
    };

  //Helper function to turn int vectors into char vectors
  auto convertTocharandRegister2 = [](NTupleReader& tr, const std::string& name)
    {
      const std::vector<bool>& boolVec = tr.getVec<bool>(name);
      std::vector<unsigned char>* charVec = new std::vector<unsigned char>(boolVec.begin(), boolVec.end());
      tr.registerDerivedVec(name+"ConvertedTochar", charVec);
      return charVec;
    };

  convertTocharandRegister1(tr, "Electron_cutBased");
  //convertTocharandRegister2(tr, "Muon_mediumId");

  //  const std::vector<float>& recoJetsBtag = *convertTofloatandRegister2(tr, "recoJetsBtag_0");
  //const std::vector<float>& qgLikelihood      = *convertTofloatandRegister2(tr, "qgLikelihood");

  //New Tagger starts here
  ttUtility::ConstAK4Inputs<float> *myConstAK4Inputs = nullptr;
  std::vector<TLorentzVector> *genTops;
  std::vector<std::vector<const TLorentzVector*>> *genTopDaughters = nullptr;
  std::vector<Constituent> constituentsMVA;
  if(tr.checkBranch("GenPart_pt"))
    {
      const std::vector<TLorentzVector>& genDecayLVec = tr.getVec_LVFromNano<float>("GenPart");
      const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("GenPart_pdgId");
      const std::vector<int>& genDecayMomIdxVec       = tr.getVec<int>("GenPart_genPartIdxMother");
      const std::vector<int>& genDecayStatFlag        = tr.getVec<int>("GenPart_statusFlags");

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
  myConstAK4Inputs->addSupplamentalVector("qgMult", *convertTofloatandRegister1(tr, "Jet_qgMult"));




  constituentsMVA = ttUtility::packageConstituents(*myConstAK4Inputs);

  //run tagger
  ttMVA->runTagger(constituentsMVA);

  delete myConstAK4Inputs;

  const TopTaggerResults& ttrMVA = ttMVA->getResults();

  const auto& candidateTops = ttrMVA.getTopCandidates();
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


  tr.registerDerivedVar("ttrMVA", &ttrMVA);

  tr.registerDerivedVar("nTops", static_cast<int>(tops.size()));

  tr.registerDerivedVec("genTops", genTops);

  tr.registerDerivedVar("highestDisc", highestDisc);

  tr.registerDerivedVar("bestTopMass", bestTopMass);
  tr.registerDerivedVar("bestTopEta", bestTopEta);
  tr.registerDerivedVar("bestTopMassLV", bestTopMassLV?(bestTopMassLV->p()):(TLorentzVector()));
  tr.registerDerivedVar("bestTopMassGenMatch", bestTopMassGenMatch);
  tr.registerDerivedVar("bestTopMassTopTag", bestTopMassTopTag);

 }
public:
PrepareTopVars(std::string taggerCfg = "TopTagger.cfg") : ttMVA(new TopTagger()) 
{

  ttMVA->setCfgFile(taggerCfg);
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

//Tagger with JEC

class PrepareTopVarsJECUp
{
 private:
  
  std::shared_ptr<TopTagger> ttMVAJECUp;
  void prepareTopVarsJECUp(NTupleReader& tr)
  {
    const std::vector<TLorentzVector>& jetsLVecTemp = tr.getVec_LVFromNano<float>("Jet");
    const std::vector<float>& recoJetsBtag      = tr.getVec<float>("Jet_btagDeepB");
    const std::vector<float>& qgLikelihood = tr.getVec<float>("Jet_qgl");  
    const std::vector<TLorentzVector> jetsLVec = tr.getVec_LVFromPtEtaPhiM<float>("Jet_pt_jesTotalUp", "Jet_eta", "Jet_phi", "Jet_mass");
    
    /*const std::vector<float>& recoJetsJecUnc = tr.getVec<float>("Jet_corr_JEC");
    for(int ijet=0; ijet<jetsLVecTemp.size(); ++ijet)
      {
	 TLorentzVector jetLVec;
	 //Up                                              
        //std:cout<<jetsLVecTemp[ijet].Pt()<<"  "<<recoJetsJecUnc[ijet]<<" M "<<jetsLVecTemp[ijet].M()<<" E "<<jetsLVecTemp[ijet].E()<<" Et "<<jetsLVecTemp[ijet].Et()<<std::endl;                                                                                                   
	 jetLVec.SetPtEtaPhiM( jetsLVecTemp[ijet].Pt()*(1 + recoJetsJecUnc[ijet]), jetsLVecTemp[ijet].Eta(), jetsLVecTemp[ijet].Phi(), jetsLVecTemp[ijet].M());
	 jetsLVec.push_back(jetLVec);
	 //std::cout<<jetLVec.Pt()<<" M "<<jetLVec.M()<<" E "<<jetLVec.E()<<" Et "<<jetLVec.Et()<<std::endl;                  
	 }*/

     auto convertTofloatandRegister1 = [](NTupleReader& tr, const std::string& name)
       {
	 const std::vector<int>& intVec = tr.getVec<int>(name);
	 std::vector<float>* floatVec = new std::vector<float>(intVec.begin(), intVec.end());
	 tr.registerDerivedVec(name+"ConvertedTofloat", floatVec);
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
	 const std::vector<int>& genDecayMomIdxVec       = tr.getVec<int>("GenPart_genPartIdxMother");
	 const std::vector<int>& genDecayStatFlag        = tr.getVec<int>("GenPart_statusFlags");
	 
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
     myConstAK4Inputs->addSupplamentalVector("qgMult", *convertTofloatandRegister1(tr, "Jet_qgMult"));
     
     
     constituentsMVA = ttUtility::packageConstituents(*myConstAK4Inputs);
     
     //run tagger                                                                                                                     
     ttMVAJECUp->runTagger(constituentsMVA);
     delete myConstAK4Inputs;
     const TopTaggerResults& ttrMVAJECUp = ttMVAJECUp->getResults();
     tr.registerDerivedVar("ttrMVAJECUp", &ttrMVAJECUp);
   }
 public:
 PrepareTopVarsJECUp(std::string taggerCfg = "TopTagger.cfg") : ttMVAJECUp(new TopTagger())
     {
       ttMVAJECUp->setCfgFile(taggerCfg);
     }
   ~PrepareTopVarsJECUp()
     {
     }
   void operator()(NTupleReader& tr)
   {
     prepareTopVarsJECUp(tr);
   }
 };


class PrepareTopVarsJECDown
{
 private:
  
  std::shared_ptr<TopTagger> ttMVAJECDown;
  void prepareTopVarsJECDown(NTupleReader& tr)
  {
    const std::vector<TLorentzVector>& jetsLVecTemp = tr.getVec_LVFromNano<float>("Jet");
    const std::vector<float>& recoJetsBtag      = tr.getVec<float>("Jet_btagDeepB");
    const std::vector<float>& qgLikelihood = tr.getVec<float>("Jet_qgl");
    const std::vector<TLorentzVector> jetsLVec = tr.getVec_LVFromPtEtaPhiM<float>("Jet_pt_jesTotalDown", "Jet_eta", "Jet_phi", "Jet_mass");
    
    /*const std::vector<float>& recoJetsJecUnc = tr.getVec<float>("Jet_corr_JEC");
    for(int ijet=0; ijet<jetsLVecTemp.size(); ++ijet)
      {
	 TLorentzVector jetLVec;
	 //Down
        //std:cout<<jetsLVecTemp[ijet].Pt()<<"  "<<recoJetsJecUnc[ijet]<<" M "<<jetsLVecTemp[ijet].M()<<" E "<<jetsLVecTemp[ijet].E()<<" Et "<<jetsLVecTemp[ijet].Et()<<std::endl;                                                                                                   
	 jetLVec.SetPtEtaPhiM( jetsLVecTemp[ijet].Pt()*(1 - recoJetsJecUnc[ijet]), jetsLVecTemp[ijet].Eta(), jetsLVecTemp[ijet].Phi(), jetsLVecTemp[ijet].M());
	 jetsLVec.push_back(jetLVec);
	 //std::cout<<jetLVec.Pt()<<" M "<<jetLVec.M()<<" E "<<jetLVec.E()<<" Et "<<jetLVec.Et()<<std::endl;                  
	 }*/

     auto convertTofloatandRegister1 = [](NTupleReader& tr, const std::string& name)
       {
	 const std::vector<int>& intVec = tr.getVec<int>(name);
	 std::vector<float>* floatVec = new std::vector<float>(intVec.begin(), intVec.end());
	 tr.registerDerivedVec(name+"ConvertedTofloat", floatVec);
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
	 const std::vector<int>& genDecayMomIdxVec       = tr.getVec<int>("GenPart_genPartIdxMother");
	 const std::vector<int>& genDecayStatFlag        = tr.getVec<int>("GenPart_statusFlags");
	 
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
     myConstAK4Inputs->addSupplamentalVector("qgMult", *convertTofloatandRegister1(tr, "Jet_qgMult"));
     
     
     constituentsMVA = ttUtility::packageConstituents(*myConstAK4Inputs);
     
     //run tagger                                                                                                                     
     ttMVAJECDown->runTagger(constituentsMVA);
     delete myConstAK4Inputs;
     const TopTaggerResults& ttrMVAJECDown = ttMVAJECDown->getResults();
     tr.registerDerivedVar("ttrMVAJECDown", &ttrMVAJECDown);
   }
 public:
 PrepareTopVarsJECDown(std::string taggerCfg = "TopTagger.cfg") : ttMVAJECDown(new TopTagger())
     {
       ttMVAJECDown->setCfgFile(taggerCfg);
     }
   ~PrepareTopVarsJECDown()
     {
     }
   void operator()(NTupleReader& tr)
   {
     prepareTopVarsJECDown(tr);
   }
 };

#endif
