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
#include "DataMC.h"

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./DataMC TTbarInc 1000 0 1" << std::endl;
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

    const unsigned int run = tr->getVar<unsigned int>("run");

    if(!isData){
      const vector<TLorentzVector> &genDecayLVec = tr->getVec_LVFromNano<float>("GenPart");
      const vector<int> &genDecayPdgIdVec = tr->getVec<int>("GenPart_pdgId");
      const vector<int> &genDecayMomIdxVec = tr->getVec<int>("GenPart_genPartIdxMother");
    }

    double ht = AnaFunctions::calcHT(jetsLVec, AnaConsts::pt20Eta24Arr);

    int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);
    int nMuons = AnaFunctions::countMuons(muonsLVec, muonsMiniIso, muonsMtw, muonsFlagIDVec, AnaConsts::muonsMiniIsoArr);
    int nElectrons = AnaFunctions::countElectrons(elesLVec, elesMiniIso, elesMtw, elesFlagIDVec, AnaConsts::elesMiniIsoArr);

    std::vector<unsigned int> *outputBIdxs = new std::vector<unsigned int>();
    const int nbJets = AnaFunctions::countCSVS(jetsLVec, recoJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr, outputBIdxs);

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
     if(isData){
       passMuTrigger = false;
       passHtTrigger = false;
       if((tr->checkBranch("HLT_IsoMu20") && tr->getVar<bool>("HLT_IsoMu20")) ||
	  (tr->checkBranch("HLT_IsoMu22") && tr->getVar<bool>("HLT_IsoMu22")) ||
	  (tr->checkBranch("HLT_IsoMu24") && tr->getVar<bool>("HLT_IsoMu24")) ||
	  (tr->checkBranch("HLT_IsoMu27") && tr->getVar<bool>("HLT_IsoMu27")) ||
	  (tr->checkBranch("HLT_IsoMu22_eta2p1") && tr->getVar<bool>("HLT_IsoMu22_eta2p1")) ||
	  (tr->checkBranch("HLT_IsoMu24_eta2p1") && tr->getVar<bool>("HLT_IsoMu24_eta2p1")) ||
	  (tr->checkBranch("HLT_IsoTkMu22") && tr->getVar<bool>("HLT_IsoTkMu22")) ||
	  (tr->checkBranch("HLT_IsoTkMu24") && tr->getVar<bool>("HLT_IsoTkMu24")) ||
	  (tr->checkBranch("HLT_Mu50")      && tr->getVar<bool>("HLT_Mu50")) ||
	  (tr->checkBranch("HLT_Mu55") && tr->getVar<bool>("HLT_Mu55"))
	  )
	 {
	   passMuTrigger = true;
	 }
	   
       if( (tr->checkBranch("HLT_PFHT750_4JetPt50") && tr->getVar<bool>("HLT_PFHT750_4JetPt50")) ||
	   (tr->checkBranch("HLT_PFHT800")          && tr->getVar<bool>("HLT_PFHT800")) ||
	   (tr->checkBranch("HLT_PFHT900")          && tr->getVar<bool>("HLT_PFHT900")) ||
	   (tr->checkBranch("HLT_PFHT1050")          && tr->getVar<bool>("HLT_PFHT1050")) ||
	   (tr->checkBranch("HLT_PFJet450") && tr->getVar<bool>("HLT_PFJet450")) ||
	   (tr->checkBranch("HLT_PFJet500") && tr->getVar<bool>("HLT_PFJet500"))
	   )
	 {
	   passHtTrigger = true;
	 }
	 
     }

     //Filters
     bool passNoiseEventFilter = true;
     bool passDataSpec = true;
     if(run >= 100000 ){ // hack to know if it's data or MC...
       auto& goodVerticesFilter = tr->getVar<bool>("Flag_goodVertices");
       // new filters
       auto& globalTightHalo2016Filter = tr->getVar<bool>("Flag_globalTightHalo2016Filter");
       auto& eeBadScFilter = tr->getVar<bool>("Flag_eeBadScFilter");
       
       passDataSpec = goodVerticesFilter && eeBadScFilter && globalTightHalo2016Filter;
     }
     
     auto& hbheNoiseFilter = tr->getVar<bool>("Flag_HBHENoiseFilter");
     auto& hbheIsoNoiseFilter = tr->getVar<bool>("Flag_HBHENoiseIsoFilter");
     auto& ecalTPFilter = tr->getVar<bool>("Flag_EcalDeadCellTriggerPrimitiveFilter");
     auto& jetID = tr->getVec<int>("Jet_jetId");
     bool jetIDFilter = true;
     for(auto& id : jetID)
       {
	 if(!(id & 0x1)) // bit 0 is for loose ID, bit 1 for tight
	   {
	     jetIDFilter = false;
	     break;
	   }
       }
     std::string type;
     tr->getType("Flag_BadPFMuonFilter", type);
     bool passBadPFMuonFilter = true;
     bool passBadChargedCandidateFilter = true;
     if(type.compare("bool") == 0)
       {
	 passBadPFMuonFilter = tr->getVar<bool>("Flag_BadPFMuonFilter");
	 passBadChargedCandidateFilter = tr->getVar<bool>("Flag_BadChargedCandidateFilter");
       }
     else if(type.compare("unsigned char") == 0)
       {
	 passBadPFMuonFilter = tr->getVar<unsigned char>("Flag_BadPFMuonFilter");
	 passBadChargedCandidateFilter = tr->getVar<unsigned char>("Flag_BadChargedCandidateFilter");
       }
     else
       {
	 passBadPFMuonFilter = tr->getVar<char>("Flag_BadPFMuonFilter");
	 passBadChargedCandidateFilter = tr->getVar<char>("Flag_BadChargedCandidateFilter");
       }
     bool passMETratioFilter = tr->getVar<float>("CaloMET_pt")!=0 ? tr->getVar<float>("MET_pt")/tr->getVar<float>("CaloMET_pt") < 5 : true;

    passNoiseEventFilter = passDataSpec && hbheNoiseFilter && hbheIsoNoiseFilter && ecalTPFilter && jetIDFilter && passBadPFMuonFilter && passBadChargedCandidateFilter && passMETratioFilter;

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

    if(!isData){
      //muontrigger
      MuTriWgt = is2016? Weight::MuTriWgt_2016(Weight::ptbin(muonL.Pt())) : is2017? Weight::MuTriWgt_2017(Weight::ptbin(muonL.Pt())) : Weight::MuTriWgt_2018(Weight::ptbin(muonL.Pt()));
      //metmhtTigger
      //HtTriWgt = Weight::HtTriWgt(Weight::metbin(met), Weight::htbin(ht));
      //pileUp weight
      PileUpWgt = tr->getVar<float>("puWeight");
      //BTag weight
      BTagWgt = tr->getVar<float>("BTagWeight");
      //TTbar weight
      bool doTTbarwgt = sample.Contains("TTbar")? true:false;
      if(doTTbarwgt){
	//TTbarWgt = tr->getVar<float>("TTbarWF");
	TTbarWgt = tr->getVar<float>("ISRWeight");
      }
      //cout<<"TTbarWgt: "<<TTbarWgt<<endl;
      tot_MCWgt_eff = MuTriWgt * PileUpWgt * BTagWgt * TTbarWgt;
      tot_MCWgt_mis = HtTriWgt * PileUpWgt * BTagWgt * TTbarWgt;
    }
    
    EventWeight = EvtWt>0? 1:-1;
    //std::cout<<"EventWeight: "<<EventWeight<<"  Lumiscale with evtwgt: "<<Lumiscale * EventWeight<<std::endl;
    //std::cout<<"mu pT: "<<muonL.Pt()<<"  is2016: "<<is2016<<"  is2017: "<<is2017<<"  is2018: "<<is2018<<" MuTriWgt: "<<MuTriWgt<<std::endl;
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
	if(recoJetsBtag[i] < 0.8) continue;//loose bjet
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

    //if((passMuonVeto!=Pass_MuonVeto)||(passEleVeto!=Pass_ElecVeto))cout<<"passMuonVeto: "<<passMuonVeto<<" Pass_MuonVeto: "<<Pass_MuonVeto<<" passEleVeto: "<<passEleVeto<<" Pass_ElecVeto: "<<Pass_ElecVeto<<endl;

    //if(passHtTrigger && passLepVeto) cout<<"passHtTrigger: "<<passHtTrigger<<" nMuons: "<<nMuons<<" nElectrons: "<<nElectrons<<" nIsoTrks: "<<nIsoTrks<<endl;     
    //cutflow
    myBaseHistgram.hcutflow_eff->Fill("Raw", 1.0);
    myBaseHistgram.hcutflow_eff->Fill("EventWeight", EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("LumiScale", Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("Filter", Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("Trigger", passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("Njet", passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("Single Muon", passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("mtW", passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("B jets", passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("#Delta #phi", passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("B Lep", passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("B Lep Mass", passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("Lep met #Delta #phi", passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("HT", passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("MET", passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("ExtraMuonVeto", passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("EleVeto", passEleVeto * passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_eff->Fill("MC trig. weights",passEleVeto * passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight * MuTriWgt);
    myBaseHistgram.hcutflow_eff->Fill("MC pileUp weights", passEleVeto * passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight * PileUpWgt);
    myBaseHistgram.hcutflow_eff->Fill("MC Btag weights", passEleVeto * passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight * BTagWgt);
    myBaseHistgram.hcutflow_eff->Fill("MC TTbar weights", passEleVeto * passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight * TTbarWgt);
    myBaseHistgram.hcutflow_eff->Fill("MC total weights", passEleVeto * passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight * tot_MCWgt_eff);
    myBaseHistgram.hcutflow_eff->Fill("HEM effect", passEleVeto * passExtraMuonVeto * passLowMET * passHT * passdPhiLep * passMassBLep * passBLep * passdPhi * passBjet * passMtW * passMu * passNjet * passMuTrigger * Pass_EventFilter * Lumiscale * EventWeight * tot_MCWgt_eff * passHEM);

    myBaseHistgram.hcutflow_mis->Fill("Raw", 1.0);
    myBaseHistgram.hcutflow_mis->Fill("EventWeight", EventWeight);
    myBaseHistgram.hcutflow_mis->Fill("LumiScale", Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_mis->Fill("Filter", Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_mis->Fill("Trigger", passHtTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_mis->Fill("Njet", passNjet * passHtTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_mis->Fill("Lepton Veto", passLepVeto * passNjet * passHtTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_mis->Fill("HT", passLepVeto * passQCDHT * passNjet * passHtTrigger * Pass_EventFilter * Lumiscale * EventWeight);
    myBaseHistgram.hcutflow_mis->Fill("MC weights", passLepVeto * passQCDHT * passNjet * passHtTrigger * Pass_EventFilter * Lumiscale * EventWeight * tot_MCWgt_mis);
    myBaseHistgram.hcutflow_mis->Fill("HEM effect", passLepVeto * passQCDHT * passNjet * passHtTrigger * Pass_EventFilter * Lumiscale * EventWeight * tot_MCWgt_mis * passHEM);

    //To compare with KenCall
    //myBaseHistgram.hcutflow_mis->Fill("Raw", 1.0);
    //myBaseHistgram.hcutflow_mis->Fill("HT trigger", passHtTrigger * 1.0);
    //myBaseHistgram.hcutflow_mis->Fill("Filter", passNoiseEventFilter * passHtTrigger * 1.0);
    //myBaseHistgram.hcutflow_mis->Fill("Lepton veto", passLepVeto * passNoiseEventFilter * passHtTrigger * 1.0);
    //myBaseHistgram.hcutflow_mis->Fill("Njet", passNjet * passLepVeto * passNoiseEventFilter * passHtTrigger * 1.0);
    //myBaseHistgram.hcutflow_mis->Fill("HT", passQCDHT * passNjet * passLepVeto * passNoiseEventFilter * passHtTrigger * 1.0);
    //myBaseHistgram.hcutflow_mis->Fill("EventWeight", passQCDHT * passNjet * passLepVeto * passNoiseEventFilter * passHtTrigger * 1.0 * EventWeight);
    //myBaseHistgram.hcutflow_mis->Fill("MC weight", passQCDHT * passNjet * passLepVeto * passNoiseEventFilter * passHtTrigger * 1.0 * EventWeight * tot_MCWgt_mis);
    //myBaseHistgram.hcutflow_mis->Fill("Lumiscale", passQCDHT * passNjet * passLepVeto * passNoiseEventFilter * passHtTrigger * 1.0 * EventWeight * tot_MCWgt_mis * Lumiscale);
    //cout<<"passNoiseEventFilter: "<<passNoiseEventFilter<<"  Pass_EventFilter: "<<Pass_EventFilter<<endl;

    if(!passNjet) continue;//to make sure resolved tagger has at least 3 jets    
    // top tagger
    //get output of tagger
    const TopTaggerResults* ttr = tr->getVar<const TopTaggerResults*>("ttrMVA");
    //Use result for top var
    vector<TopObject*> Ntop = ttr->getTops();
    vector<TopObject> NtopCand = ttr->getTopCandidates();
    //    std::cout<<"Ntop: "<<Ntop.size()<<" NtopCand: "<<NtopCand.size()<<std::endl;
    pUtility::FillInt(myBaseHistgram.hNtop,Ntop.size(),Lumiscale * EventWeight);
    pUtility::FillInt(myBaseHistgram.hNtopCand,NtopCand.size(),Lumiscale * EventWeight);

    for(const TopObject& topcand : NtopCand)
      {
	//std::cout<<"Cand disc: "<<topcand.getDiscriminator()<<std::endl;
	pUtility::FillDouble(myBaseHistgram.htopCandDisc,topcand.getDiscriminator(),Lumiscale * EventWeight);
      }


    //Efficiency & Mistag Rate are from same formula but derived with different CS:
    const double discThrLWP = 0.75; 
    const double discThrMWP = 0.85; 
    const double discThrTWP = 0.95; 
    const double discThrXWP = 0.92;
    TopObject bestTopCand = GetBestTopCand(NtopCand);
    bool isbestCandTag = CandTag(bestTopCand, Ntop);
    int genmatch = isData? 0:std::min(partonMatch(bestTopCand, 0.6), constMatch(bestTopCand, 0.6));

    //For SF validation
    bool isSFVal = isData? false: true;
    float SF = 1.0;
    if(isSFVal){
      for(const auto& top : Ntop)
	{
	  TopObject Top = *top;
	  int match = std::min(partonMatch(Top, 0.6), constMatch(Top, 0.6));
	  float sf = (match==3 || match==2)? Weight::TopEffSF_17 : Weight::TopMisSF_17;
	  //float sf = (match==3 || match==2)? Weight::TopEffSF_2016(Weight::Topptbin(Top.p().Pt()), Weight::jetbin(cntNJetsPt30Eta24)) : Weight::TopMisSF_2016(Weight::Topptbin(Top.p().Pt()), Weight::jetbin(cntNJetsPt30Eta24));
	  //cout<<"match: "<<match<<" sf: "<<sf<<endl;
	  SF = SF*sf;
	}
    }
    //cout<<" SF: "<<SF<<endl;
    float tot_MCWgt_effSF = 1.0;
    float tot_MCWgt_misSF = 1.0;
    tot_MCWgt_effSF = tot_MCWgt_eff * SF;
    tot_MCWgt_misSF = tot_MCWgt_mis * SF;
    

    //Eff
    //Single muon CS (data triggered with single mu)
    if(passNjet && passMu && passBjet && passdPhi && passBLep && passMassBLep && passdPhiLep && passMtW && passHT && passLowMET && passExtraMuonVeto && passEleVeto && passMuTrigger && Pass_EventFilter && passHEM){
    //if(passNjet && passMu && passBjet && passdPhi && passdPhiLep && passMtW && passHT && passLowMET && passMuTrigger && passNoiseEventFilter){
      //Some distribution with muon trigger weight
      pUtility::FillInt(myBaseHistgram.hNjet_Eff, cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hmet_Eff, met, Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hht_Eff, ht, Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hBestCandDisc_Eff, bestTopCand.getDiscriminator(),Lumiscale * EventWeight*tot_MCWgt_eff);     
      myBaseHistgram.hNjetDisc_Eff->Fill(cntNJetsPt30Eta24, bestTopCand.getDiscriminator(), Lumiscale * EventWeight*tot_MCWgt_eff);

      //den
      //if(bestTopCand.p().M()==0.0) cout<<"bestTopCand.p().Pt(): "<<bestTopCand.p().Pt()<<"  bestTopCand.getDiscriminator(): "<<bestTopCand.getDiscriminator()<<endl;
      if(bestTopCand.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	pUtility::Fill2D(myBaseHistgram.htopPtNjetEff_den,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_eff);
	//num
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetEff_num_LWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_effSF);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetEff_num_MWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_effSF);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetEff_num_TWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_effSF);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetEff_num_XWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_effSF);
	}
	//Extra plot with gen matching #for MC
	pUtility::FillDouble(myBaseHistgram.hrecotopPt,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_0m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_1m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_2m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_3m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
      
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_0m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_1m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_2m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_3m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_0m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_1m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_2m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_3m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_0m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_1m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_2m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_3m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_0m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_1m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_2m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_3m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}

      }
      //for events variables
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hmetEff_den,met,Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hhtEff_den,ht,Lumiscale * EventWeight*tot_MCWgt_eff);
      //num
      if(isTopPassDisc(Ntop,discThrLWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hmetEff_num_LWP,met,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hhtEff_num_LWP,ht,Lumiscale * EventWeight*tot_MCWgt_effSF);
	
	}
      if(isTopPassDisc(Ntop,discThrMWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hmetEff_num_MWP,met,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hhtEff_num_MWP,ht,Lumiscale * EventWeight*tot_MCWgt_effSF);
	}
      if(isTopPassDisc(Ntop,discThrTWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hmetEff_num_TWP,met,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hhtEff_num_TWP,ht,Lumiscale * EventWeight*tot_MCWgt_effSF);
	}
      if(isTopPassDisc(Ntop,discThrXWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hmetEff_num_XWP,met,Lumiscale * EventWeight*tot_MCWgt_effSF);
	  pUtility::FillDouble(myBaseHistgram.hhtEff_num_XWP,ht,Lumiscale * EventWeight*tot_MCWgt_effSF);
	}
      //Extra plot with gen matching #for MC 
      if(bestTopCand.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.hNjetEff,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	if(!isData && (partonMatch(bestTopCand, 0.6)==-1 || partonMatch(bestTopCand, 0.6)==0)) pUtility::FillDouble(myBaseHistgram.hNjetEff_0m,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	if(!isData && (partonMatch(bestTopCand, 0.6)==1)) pUtility::FillDouble(myBaseHistgram.hNjetEff_1m,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff); 
	if(!isData && (partonMatch(bestTopCand, 0.6)==2)) pUtility::FillDouble(myBaseHistgram.hNjetEff_2m,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff); 
	if(!isData && (partonMatch(bestTopCand, 0.6)==3)) pUtility::FillDouble(myBaseHistgram.hNjetEff_3m,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff); 
	
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==-1 || partonMatch(bestTopCand, 0.6)==0)) pUtility::FillDouble(myBaseHistgram.hNjetEff_0m_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==1)) pUtility::FillDouble(myBaseHistgram.hNjetEff_1m_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==2)) pUtility::FillDouble(myBaseHistgram.hNjetEff_2m_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==3)) pUtility::FillDouble(myBaseHistgram.hNjetEff_3m_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==-1 || partonMatch(bestTopCand, 0.6)==0)) pUtility::FillDouble(myBaseHistgram.hNjetEff_0m_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==1)) pUtility::FillDouble(myBaseHistgram.hNjetEff_1m_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==2)) pUtility::FillDouble(myBaseHistgram.hNjetEff_2m_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==3)) pUtility::FillDouble(myBaseHistgram.hNjetEff_3m_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==-1 || partonMatch(bestTopCand, 0.6)==0)) pUtility::FillDouble(myBaseHistgram.hNjetEff_0m_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==1)) pUtility::FillDouble(myBaseHistgram.hNjetEff_1m_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==2)) pUtility::FillDouble(myBaseHistgram.hNjetEff_2m_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==3)) pUtility::FillDouble(myBaseHistgram.hNjetEff_3m_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==-1 || partonMatch(bestTopCand, 0.6)==0)) pUtility::FillDouble(myBaseHistgram.hNjetEff_0m_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==1)) pUtility::FillDouble(myBaseHistgram.hNjetEff_1m_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==2)) pUtility::FillDouble(myBaseHistgram.hNjetEff_2m_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && (partonMatch(bestTopCand, 0.6)==3)) pUtility::FillDouble(myBaseHistgram.hNjetEff_3m_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      }
      
    }

    

    //Mistag Rate
    //QCD CS (data triggered with with ht)
    if(passNjet && passLepVeto && passQCDHT && passHtTrigger && Pass_EventFilter && passHEM){
      //Some distribution with ht trigger weight
      pUtility::FillInt(myBaseHistgram.hNjet_Mis, cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_mis);
      pUtility::FillDouble(myBaseHistgram.hmet_Mis, met, Lumiscale * EventWeight*tot_MCWgt_mis);
      pUtility::FillDouble(myBaseHistgram.hht_Mis, ht, Lumiscale * EventWeight*tot_MCWgt_mis);
      pUtility::FillDouble(myBaseHistgram.hBestCandDisc_Mis, bestTopCand.getDiscriminator(),Lumiscale * EventWeight*tot_MCWgt_mis);     
      myBaseHistgram.hNjetDisc_Mis->Fill(cntNJetsPt30Eta24, bestTopCand.getDiscriminator(), Lumiscale * EventWeight*tot_MCWgt_mis);
      //for best cand.                                                                                                                        
      //den                                                                                                                                   
      if(bestTopCand.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	pUtility::Fill2D(myBaseHistgram.htopPtNjetMis_den,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_mis);
        //num                                                                                                                                 
        if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetMis_num_LWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_misSF);
        }
        if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetMis_num_MWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_misSF);
        }
        if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetMis_num_TWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_misSF);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::Fill2D(myBaseHistgram.htopPtNjetMis_num_XWP,bestTopCand.p().Pt(), cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_misSF);
        }
      	//Extra plot with gen matching #for MC
	pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_0m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_1m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_2m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_3m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
      
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_0m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_1m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_2m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_3m_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_0m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_1m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_2m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_3m_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_0m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_1m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_2m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_3m_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && (genmatch==-1 || genmatch==0)) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_0m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff); 
	  if(!isData && genmatch==1) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_1m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==2) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_2m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  if(!isData && genmatch==3) pUtility::FillDouble(myBaseHistgram.hrecotopPt_Mis_3m_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      }
      //for events variables
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      pUtility::FillDouble(myBaseHistgram.hmetMis_den,met,Lumiscale * EventWeight*tot_MCWgt_mis);
      pUtility::FillDouble(myBaseHistgram.hhtMis_den,ht,Lumiscale * EventWeight*tot_MCWgt_mis);

      //num
      if(isTopPassDisc(Ntop,discThrLWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hmetMis_num_LWP,met,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hhtMis_num_LWP,ht,Lumiscale * EventWeight*tot_MCWgt_misSF);
	}
      if(isTopPassDisc(Ntop,discThrMWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hmetMis_num_MWP,met,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hhtMis_num_MWP,ht,Lumiscale * EventWeight*tot_MCWgt_misSF);	
	}
      if(isTopPassDisc(Ntop,discThrTWP))
	{
       	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hmetMis_num_TWP,met,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hhtMis_num_TWP,ht,Lumiscale * EventWeight*tot_MCWgt_misSF);
	}
      if(isTopPassDisc(Ntop,discThrXWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hmetMis_num_XWP,met,Lumiscale * EventWeight*tot_MCWgt_misSF);
	  pUtility::FillDouble(myBaseHistgram.hhtMis_num_XWP,ht,Lumiscale * EventWeight*tot_MCWgt_misSF);
	}

    }
    

    //Dilep Control region
    if(passDilep && passMedMET && passHT && passMuTrigger && Pass_EventFilter && passHEM){
      pUtility::FillInt(myBaseHistgram.hNjet_2l, cntNJetsPt30Eta24, Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hmet_2l, met, Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hht_2l, ht, Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hBestCandDisc_2l, bestTopCand.getDiscriminator(),Lumiscale * EventWeight*tot_MCWgt_eff);     
      if(bestTopCand.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPt2l_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	//num
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPt2l_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPt2l_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPt2l_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPt2l_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      }
      //for events variables
      //den
      pUtility::FillDouble(myBaseHistgram.hNjet2l_den,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
      //num
      if(isTopPassDisc(Ntop,discThrLWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjet2l_num_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(Ntop,discThrMWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjet2l_num_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(Ntop,discThrTWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjet2l_num_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(Ntop,discThrXWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjet2l_num_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}

    }

  }//event loop
  (myBaseHistgram.oFile)->Write();
  return 0;
}
