#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/PileupWeights.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/TTbarCorrector.h"
#include "SusyAnaTools/Tools/ISRCorrector.h"

#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"

#include "ResTagger.h"
#include "PlotUtility.h"
#include "Weight.h"
#include "Eff.h"

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Eff TTbarInc 1000 0 1" << std::endl;
      return -1;
    }
  const char *subsamplename = argv[1]; 
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);  
  const int maxevent = std::atoi(Maxevent);
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  const string condorSpec = argc==6 ? argv[5]: "";  
  bool iscondor = condorSpec.empty()? false:true;
  TChain *fChain = 0;
  if(!FillChain(fChain, subsamplename, iscondor, startfile, filerun))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }
  

  std::string spec = "TagVal";
  TString sample(subsamplename); 
  bool isData(false);
  NTupleReader *tr = nullptr;
  tr = new NTupleReader(fChain, {"run"});

  //top tagger  
  PrepareTopVars prepareTopVars("TopTagger.cfg");
  tr->registerFunction(prepareTopVars);

  TopCat topcat;
  TopVar topvar;
  if(sample.Contains("Data") || sample.Contains("DATA")) isData = true;
  bool is2016 = sample.Contains("2016")? true:false;
  bool is2017 = sample.Contains("2017")? true:false;
  bool is2018 = sample.Contains("2018")? true:false;

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;
  int ev = 0;
  // Loop over the events (tree entries)
  while(tr->getNextEvent()){
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;
    ev++;

    const vector<TLorentzVector> &muonsLVec = tr->getVec_LVFromNano<float>("Muon");
    const vector<float> &muonsMiniIso = tr->getVec<float>("Muon_miniPFRelIso_all");
    const vector<unsigned char> & muonsFlagIDVec = tr->getVec<unsigned char>("Muon_mediumId");
    const vector<int> & muonsChrg = tr->getVec<int>("Muon_charge");
    
    const vector<TLorentzVector> &elesLVec = tr->getVec_LVFromNano<float>("Electron");
    const vector<float> &elesMiniIso = tr->getVec<float>("Electron_miniPFRelIso_all");
    const vector<unsigned char> & elesFlagIDVec = tr->getVec<unsigned char>("Electron_cutBasedConvertedTochar");
    const vector<int> & elesChrg = tr->getVec<int>("Electron_charge");


    const vector<TLorentzVector> &jetsLVec = tr->getVec_LVFromNano<float>("Jet");
    const vector<float> &recoJetsBtag = tr->getVec<float>("Jet_btagDeepB");
    
    
    float met=tr->getVar<float>("MET_pt");
    float metphi=tr->getVar<float>("MET_phi");
    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

    vector<float> muonsMtw(muonsLVec.size());
    for(int i = 0; i < muonsMtw.size(); ++i)
      {
	muonsMtw[i] = sqrt(2*metLVec.Pt()*muonsLVec[i].Pt()*(1-cos(ROOT::Math::VectorUtil::DeltaPhi(metLVec, muonsLVec[i]))));
      }
    vector<float> elesMtw(elesLVec.size());
    for(int i = 0; i < elesMtw.size(); ++i)
      {
	elesMtw[i] = sqrt(2*metLVec.Pt()*elesLVec[i].Pt()*(1-cos(ROOT::Math::VectorUtil::DeltaPhi(metLVec, elesLVec[i]))));
      }

    float EvtWt = isData? 1 : tr->getVar<float>("Generator_weight");
    //std::cout<<"EvtWt"<<EvtWt<<std::endl;
    EventWeight = EvtWt>0? 1:-1;
    const unsigned int run = tr->getVar<unsigned int>("run");
    if(!isData){
      const vector<TLorentzVector> &genDecayLVec = tr->getVec_LVFromNano<float>("GenPart");
      const vector<int> &genDecayPdgIdVec = tr->getVec<int>("GenPart_pdgId");
      const vector<int> &genDecayMomIdxVec = tr->getVec<int>("GenPart_genPartIdxMother");
      const vector<int>& genDecayStatFlag        = tr->getVec<int>("GenPart_statusFlags");
   
      //gen plots for signal 
      auto gentopinfo = ttUtility::GetTopdauGenLVecFromNano(genDecayLVec, genDecayPdgIdVec, genDecayStatFlag, genDecayMomIdxVec);
      vector<TLorentzVector> gentops = gentopinfo.first;
      for(const auto& top:gentops){
	pUtility::FillDouble(myBaseHistgram.hgentopPt, top.Pt(), Lumiscale * EventWeight);
      }
      bool neu100 = false;
      for(unsigned ig=0; ig < genDecayLVec.size(); ig++)
	{
	  int pdgId = genDecayPdgIdVec[ig];
	  if(pdgId == 1000006) pUtility::FillDouble(myBaseHistgram.hstopMass, genDecayLVec[ig].M(), Lumiscale * EventWeight);
	  if(pdgId > 1000021) pUtility::FillDouble(myBaseHistgram.hneutralinoMass, genDecayLVec[ig].M(), Lumiscale * EventWeight);
	  if(pdgId > 1000021 && genDecayLVec[ig].M()==100) neu100 = true;
	}
    }

    double ht = AnaFunctions::calcHT(jetsLVec, AnaConsts::pt20Eta24Arr);

    int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);
    int nMuons = AnaFunctions::countMuons(muonsLVec, muonsMiniIso, muonsMtw, muonsFlagIDVec, AnaConsts::muonsMiniIsoArr);
    int nElectrons = AnaFunctions::countElectrons(elesLVec, elesMiniIso, elesMtw, elesFlagIDVec, AnaConsts::elesMiniIsoArr);

    string key = is2016? "2016MC" : is2017? "2017MC" : "2018MC";
    string bwpM = "cutM";
    string bwpL = "cutL";

    map<string, map<string, float>> DS = AnaConsts::DeepCSV;
    map<string, float> DC = DS[key];
    float cutDeepCSV = DC[bwpM];
    float looseDeepCSV = DC[bwpL];
    //cout<<"key: "<<key.c_str()<<"  cutDeepCSV: "<<cutDeepCSV<<endl;

    std::vector<unsigned int> *outputBIdxs = new std::vector<unsigned int>();
    const int nbJets = AnaFunctions::countCSVS(jetsLVec, recoJetsBtag, cutDeepCSV, AnaConsts::bTagArr, outputBIdxs);

    int nIsoTrks = 0;

    //dPhi calc
    std::vector<float> deltaPhiVec = AnaFunctions::calcDPhi(jetsLVec, metLVec, 3, AnaConsts::dphiArr);    
    //Get Lepton
    float MupT = 40;
    std::vector<int> passMuIdx = GetMuIdx(muonsLVec, muonsMiniIso, muonsMtw, muonsFlagIDVec, AnaConsts::muonsMiniIsoArr);
    int selIdx(-1);
    TLorentzVector muonL = GetMuLVec(muonsLVec, passMuIdx, selIdx, MupT);
    //cout<<"muonsLVec: "<<muonsLVec.size()<<"  passMuIdx: "<<passMuIdx.size()<<"Mu pT: "<<muonL.Pt()<<"  selid: "<<selIdx<<endl;
    std::vector<TLorentzVector> ExtraMuLVec = GetLep(muonsLVec, passMuIdx, selIdx);
    //cout<<"ExtraMuLVec: "<<ExtraMuLVec.size()<<endl;    

    //Trigger
    bool passMuTrigger = true;
    bool passHtTrigger = true; 

    
    //boolian from ntuple
    const bool Pass_EventFilter = tr->getVar<bool>("Pass_EventFilter");
    const bool Pass_ElecVeto = tr->getVar<bool>("Pass_ElecVeto");
    const bool Pass_MuonVeto = tr->getVar<bool>("Pass_MuonVeto");
    const bool Pass_IsoTrkVeto = tr->getVar<bool>("Pass_IsoTrkVeto");

    //MC weights
    float HtTriWgt = 1.0;
    float MuTriWgt = 1.0;
    float BTagWgt = 1.0;
    float PileUpWgt = 1.0;
    float TTbarWgt = 1.0;

    float tot_MCWgt_eff = 1.0;
    float tot_MCWgt_mis = 1.0;

    
    

    //cut
    bool passNjet(false);
    bool passMu(false);
    bool passBjet(false);
    bool passdPhi(false);
    bool passBLep(false);
    bool passMassBLep(false);
    bool passdPhiLep(false);
    bool passMtW(false);
    bool passHT(false);
    bool passHighMET(false);
    bool passLowMET(false);
    bool passExtraMuonVeto(false);
    bool passEleVeto(false);
    bool passMuonVeto(false);
    bool passQCDHT(false);
    bool passLepVeto(false);
    bool passMedMET(false);
    bool passDilep = CheckDiLep(muonsLVec, muonsMiniIso, muonsMtw, muonsFlagIDVec, AnaConsts::muonsMiniIsoArr, elesLVec, elesMiniIso, elesMtw, elesFlagIDVec, AnaConsts::elesMiniIsoArr, muonsChrg, elesChrg, MupT, 20);
    bool passHEM = is2018? PassHEM(jetsLVec, -3.2, -1.2, -1.77, -0.67, 30): true;
    if(cntNJetsPt30Eta24>=4) passNjet = true;
    if(muonL.Pt()) passMu = true;
    if(selIdx>=0){if(muonsMtw[selIdx]<100)passMtW = true;}
    if(nbJets) passBjet = true;
    if( deltaPhiVec.at(0) > AnaConsts::dPhi0_CUT && deltaPhiVec.at(1) > AnaConsts::dPhi1_CUT && deltaPhiVec.at(2) > AnaConsts::dPhi2_CUT){
      passdPhi = true;
    }

    for(int i = 0; i < jetsLVec.size(); i++)
      {
	if(recoJetsBtag[i] < looseDeepCSV) continue;//loose bjet
	if(jetsLVec[i].DeltaR(muonL) < 1.5)
	  {
	    passBLep = true;
	  }
	if((muonL + jetsLVec[i]).M()>30 && (muonL + jetsLVec[i]).M()<180)
	  {
	    passMassBLep = true;
	  }
      }
    
    if(fabs(muonL.DeltaPhi(metLVec))<0.8)passdPhiLep = true;
    if(ht>200)passHT = true;
    if(met>50)passLowMET = true;
    if(met>100)passMedMET = true;
    if(met>250)passHighMET = true;
    if(ExtraMuLVec.size()==0)passExtraMuonVeto = true;
    if(nElectrons==0)passEleVeto = true;
    if(nMuons==0)passMuonVeto = true;

    if(nMuons==0 && nElectrons==0 && Pass_IsoTrkVeto) passLepVeto = true;
    float HTcut = is2016?1000:1200;
    if(ht>HTcut)passQCDHT=true;


    if(!passNjet) continue;//to make sure resolved tagger has at least 3 jets    
    // top tagger
    //get output of tagger
    const TopTaggerResults* ttr = tr->getVar<const TopTaggerResults*>("ttrMVA");
    //Use result for top var
    vector<TopObject*> Ntop = ttr->getTops();
    vector<TopObject> NtopCand = ttr->getTopCandidates();
    //    std::cout<<"Ntop: "<<Ntop.size()<<" NtopCand: "<<NtopCand.size()<<std::endl;

    //Efficiency & Mistag Rate are from same formula but derived with different CS:
    const double discThrLWP = 0.75; 
    const double discThrMWP = 0.85; 
    const double discThrTWP = 0.95; 
    const double discThrXWP = 0.92;
    TopObject bestTopCand = GetBestTopCand(NtopCand);
    bool isbestCandTag = CandTag(bestTopCand, Ntop);
    bool isconstB = CheckBconst(bestTopCand, cutDeepCSV);


    //Eff
    //Single muon CS (data triggered with single mu)
    //if(passNjet && passMu && passBjet && passdPhi && passBLep && passMassBLep && passdPhiLep && passMtW && passHT && passLowMET && passExtraMuonVeto && passEleVeto && passMuTrigger && Pass_EventFilter && passHEM){
    if(!isData){
      int genmatch = std::min(partonMatch(bestTopCand, 0.6), constMatch(bestTopCand, 0.6));    
      if(passNjet && passLowMET && Pass_EventFilter){
	if(bestTopCand.p().M()>0.0){
	  //den
	  if(genmatch==2 || genmatch==3) {
	    pUtility::FillDouble(myBaseHistgram.hrecotopPt,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	      pUtility::FillDouble(myBaseHistgram.hrecotopPt_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	    }
	    if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	      pUtility::FillDouble(myBaseHistgram.hrecotopPt_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	    }
	    if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	      pUtility::FillDouble(myBaseHistgram.hrecotopPt_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	    }
	    if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	      pUtility::FillDouble(myBaseHistgram.hrecotopPt_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	    }
	  }
	}
      }//CS
    }
    //QCD CS (mistag rate check in presence of b)
    if(passNjet && passLepVeto && passQCDHT && passHtTrigger && Pass_EventFilter && passHEM){
      //Top pT
      if(bestTopCand.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtmis_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	if(isconstB)pUtility::FillDouble(myBaseHistgram.htopPtmisB_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	else pUtility::FillDouble(myBaseHistgram.htopPtmisNoB_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtmis_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  if(isconstB)pUtility::FillDouble(myBaseHistgram.htopPtmisB_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  else pUtility::FillDouble(myBaseHistgram.htopPtmisNoB_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtmis_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  if(isconstB)pUtility::FillDouble(myBaseHistgram.htopPtmisB_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  else pUtility::FillDouble(myBaseHistgram.htopPtmisNoB_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtmis_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  if(isconstB)pUtility::FillDouble(myBaseHistgram.htopPtmisB_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  else pUtility::FillDouble(myBaseHistgram.htopPtmisNoB_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtmis_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  if(isconstB)pUtility::FillDouble(myBaseHistgram.htopPtmisB_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	  else pUtility::FillDouble(myBaseHistgram.htopPtmisNoB_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis); 
	}
      }
      //Njet
      pUtility::FillDouble(myBaseHistgram.hNjetmis_den,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      if(isTopPassDisc(Ntop,discThrLWP)){
	pUtility::FillDouble(myBaseHistgram.hNjetmis_num_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      }
      if(isTopPassDisc(Ntop,discThrMWP)){
	pUtility::FillDouble(myBaseHistgram.hNjetmis_num_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      }
      if(isTopPassDisc(Ntop,discThrTWP)){
	pUtility::FillDouble(myBaseHistgram.hNjetmis_num_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      }
      if(isTopPassDisc(Ntop,discThrXWP)){
	pUtility::FillDouble(myBaseHistgram.hNjetmis_num_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      }

    }
    //Dilep Control region
    if(passDilep && passMedMET && passHT && passMuTrigger && Pass_EventFilter && passHEM){
    }


  }//event loop
  (myBaseHistgram.oFile)->Write();
  return 0;
}
