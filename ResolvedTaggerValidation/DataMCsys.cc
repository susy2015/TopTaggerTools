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
#include "DataMCsys.h"

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./DataMCsys TTbarInc 1000 0 1" << std::endl;
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
  bool isGenSys(false);
  NTupleReader *tr =0;
  tr = new NTupleReader(fChain, {"run"});

  if(sample.Contains("Data") || sample.Contains("DATA")) isData = true;
  if(sample.Contains("Powheg") || sample.Contains("Herwigpp") || sample.Contains("_herwigpp") || sample.Contains("_pythia")) isGenSys = true;
  bool is2016 = sample.Contains("2016")? true:false;
  bool is2017 = sample.Contains("2017")? true:false;
  bool is2018 = sample.Contains("2018")? true:false;

  // cout<<"isGenSys: "<<isGenSys<<endl;
  //top tagger  
  PrepareTopVars prepareTopVars("TopTagger.cfg");
  tr->registerFunction(prepareTopVars);
  //JEC
  if(!isGenSys){
    PrepareTopVarsJECUp prepareTopVarsJECUp("TopTagger.cfg");
    tr->registerFunction(prepareTopVarsJECUp);
    PrepareTopVarsJECDown prepareTopVarsJECDown("TopTagger.cfg");
    tr->registerFunction(prepareTopVarsJECDown);
  }

  TopCat topcat;
  TopVar topvar;


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
    // std::cout<<"Event: "<<ev<<std::endl;

    const vector<TLorentzVector> &muonsLVec = tr->getVec_LVFromNano<float>("Muon");
    const vector<float> &muonsMiniIso = tr->getVec<float>("Muon_miniPFRelIso_all");
    const vector<unsigned char> & muonsFlagIDVec = tr->getVec<unsigned char>("Muon_mediumId");
   
    const vector<TLorentzVector> &elesLVec = tr->getVec_LVFromNano<float>("Electron");
    const vector<float> &elesMiniIso = tr->getVec<float>("Electron_miniPFRelIso_all");
    const vector<unsigned char> & elesFlagIDVec = tr->getVec<unsigned char>("Electron_cutBasedConvertedTochar");
   
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
    //std::cout<<"EvtWt: "<<EvtWt<<std::endl;
    //double EvtWt = 1.0;
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
    int nIsoTrks = 0;
    std::vector<unsigned int> *outputBIdxs = new std::vector<unsigned int>();
    const int nbJets = AnaFunctions::countCSVS(jetsLVec, recoJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr, outputBIdxs);

    
    //dPhi calc
    std::vector<float> deltaPhiVec = AnaFunctions::calcDPhi(jetsLVec, metLVec, 3, AnaConsts::dphiArr);    
    //Get Lepton
    std::vector<int> passMuIdx = GetMuIdx(muonsLVec, muonsMiniIso, muonsMtw, muonsFlagIDVec, AnaConsts::muonsMiniIsoArr);
    int selIdx(-1);
    TLorentzVector muonL = GetMuLVec(muonsLVec, passMuIdx, selIdx, 40.0);
    //cout<<"muonsLVec: "<<muonsLVec.size()<<"  passMuIdx: "<<passMuIdx.size()<<"Mu pT: "<<muonL.Pt()<<"  selid: "<<selIdx<<endl;
    std::vector<TLorentzVector> ExtraMuLVec = GetLep(muonsLVec, passMuIdx, selIdx);
    //cout<<"ExtraMuLVec: "<<ExtraMuLVec.size()<<endl;
    
    //Trigger
     bool passMuTrigger = true;
     bool passHtTrigger = true; 
     if(isData){
       passMuTrigger = false;
       passHtTrigger = false;
       if((tr->checkBranch("HLT_IsoMu24") && tr->getVar<bool>("HLT_IsoMu24")) ||
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
           (tr->checkBranch("HLT_PFJet450") && tr->getVar<bool>("HLT_PFJet450"))
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


    //Sys
    float BTagWgtUp = 1.0;
    float BTagWgtDown = 1.0;
    float PileUpWgtUp = 1.0;
    float PileUpWgtDown = 1.0;
    float PDFWgtUp = 1.0;
    float PDFWgtDown = 1.0;

    float tot_MCWgt_eff_BTUp = 1.0;
    float tot_MCWgt_eff_BTDown = 1.0;
    float tot_MCWgt_mis_BTUp = 1.0;
    float tot_MCWgt_mis_BTDown = 1.0;
    float tot_MCWgt_eff_PUUp = 1.0;
    float tot_MCWgt_eff_PUDown = 1.0;
    float tot_MCWgt_mis_PUUp = 1.0;
    float tot_MCWgt_mis_PUDown = 1.0;
    float tot_MCWgt_eff_PDFUp = 1.0;
    float tot_MCWgt_eff_PDFDown = 1.0;
    float tot_MCWgt_mis_PDFUp = 1.0;
    float tot_MCWgt_mis_PDFDown = 1.0;



    if(!isData){
      //muontrigger
	MuTriWgt = is2016? Weight::MuTriWgt_2016(Weight::ptbin(muonL.Pt())) : is2017? Weight::MuTriWgt_2017(Weight::ptbin(muonL.Pt())) : Weight::MuTriWgt_2018(Weight::ptbin(muonL.Pt()));
      //metmhtTigger
      //HtTriWgt = Weight::HtTriWgt(Weight::metbin(met), Weight::htbin(ht));
      //pileUp weight
      PileUpWgt = tr->getVar<float>("puWeight");
      PileUpWgtUp = tr->getVar<float>("puWeight_Up");
      PileUpWgtDown = tr->getVar<float>("puWeight_Down");
      //BTag weight
      BTagWgt = tr->getVar<float>("BTagWeight");
      BTagWgtUp = tr->getVar<float>("BTagWeight_Up");
      BTagWgtDown = tr->getVar<float>("BTagWeight_Down");
      
      //TTbar weight
      bool doTTbarwgt = sample.Contains("TTbar")? true:false;
      if(doTTbarwgt && !isGenSys){
	//TTbarWgt = tr->getVar<float>("TTbarWF");
	TTbarWgt = tr->getVar<float>("ISRWeight");
      }
      
      //std::cout<<"TTbarWgt: "<<TTbarWgt<<std::endl;
      //PDF Unc weight
      PDFWgtUp = tr->getVar<float>("pdfWeight_Up");
      PDFWgtDown = tr->getVar<float>("pdfWeight_Down");
      
      //Total
      tot_MCWgt_eff = MuTriWgt * PileUpWgt * BTagWgt * TTbarWgt;
      tot_MCWgt_mis = HtTriWgt * PileUpWgt * BTagWgt * TTbarWgt;
      
      tot_MCWgt_eff_BTUp = MuTriWgt * PileUpWgt * BTagWgtUp * TTbarWgt;
      tot_MCWgt_mis_BTUp = HtTriWgt * PileUpWgt * BTagWgtUp * TTbarWgt;
      tot_MCWgt_eff_BTDown = MuTriWgt * PileUpWgt * BTagWgtDown * TTbarWgt;
      tot_MCWgt_mis_BTDown = HtTriWgt * PileUpWgt * BTagWgtDown * TTbarWgt;

      tot_MCWgt_eff_PUUp = MuTriWgt * PileUpWgtUp * BTagWgt * TTbarWgt;
      tot_MCWgt_mis_PUUp = HtTriWgt * PileUpWgtUp * BTagWgt * TTbarWgt;
      tot_MCWgt_eff_PUDown = MuTriWgt * PileUpWgtDown * BTagWgt * TTbarWgt;
      tot_MCWgt_mis_PUDown = HtTriWgt * PileUpWgtDown * BTagWgt * TTbarWgt;

      tot_MCWgt_eff_PDFUp = MuTriWgt * PileUpWgt * BTagWgt * TTbarWgt * PDFWgtUp;
      tot_MCWgt_mis_PDFUp = HtTriWgt * PileUpWgt * BTagWgt * TTbarWgt * PDFWgtUp;
      tot_MCWgt_eff_PDFDown = MuTriWgt * PileUpWgt * BTagWgt * TTbarWgt * PDFWgtDown;
      tot_MCWgt_mis_PDFDown = HtTriWgt * PileUpWgt * BTagWgt * TTbarWgt * PDFWgtDown;

    }
    if(isGenSys) EventWeight = EvtWt;
    else EventWeight = EvtWt>0? 1:-1;    
    //std::cout<<"EventWeight: "<<EventWeight<<std::endl;
    //std::cout<<"EventWeight: "<<EventWeight<<"  Lumiscale with evtwgt: "<<Lumiscale * EventWeight<<std::endl;
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

    bool passLepVeto(false);
    bool passQCDHT(false);

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
    if(met>250)passHighMET = true;
    if(ExtraMuLVec.size()==0)passExtraMuonVeto = true;
    if(nElectrons==0)passEleVeto = true;

    if(nMuons==0 && nElectrons==0 && Pass_IsoTrkVeto==0) passLepVeto = true;
    float HTcut = is2016?1000:1200;
    if(ht>HTcut)passQCDHT=true;


    if(!passNjet) continue;//to make sure resolved tagger has at least 3 jets    
    // top tagger
    //get output of tagger
    const TopTaggerResults* ttr = tr->getVar<const TopTaggerResults*>("ttrMVA");
    //Use result for top var
    vector<TopObject*> Ntop = ttr->getTops();
    vector<TopObject> NtopCand = ttr->getTopCandidates();
    //std::cout<<"Nominal: "<<"Ntop: "<<Ntop.size()<<" NtopCand: "<<NtopCand.size()<<std::endl;

    //Efficiency & Mistag Rate are from same formula but derived with different CS:
    const double discThrLWP = 0.75; 
    const double discThrMWP = 0.85; 
    const double discThrTWP = 0.95; 
    const double discThrXWP = 0.92;
    TopObject bestTopCand = GetBestTopCand(NtopCand);
    bool isbestCandTag = CandTag(bestTopCand, Ntop);

    //Eff
    //Single muon CS (data reiggered with single mu)
    //cout<<"passNjet: "<<passNjet<<" passMu: "<<passMu<<" passBjet  "<<passBjet<<" passdPhi  "<<passdPhi<<" passBLep "<<passBLep<<" passdPhiLep  "<<passdPhiLep<<" Pass_EventFilter  "<<Pass_EventFilter<<endl;
    if(passNjet && passMu && passBjet && passdPhi && passBLep && passMassBLep && passdPhiLep && passMtW && passHT && passLowMET && passExtraMuonVeto && passEleVeto && passMuTrigger && Pass_EventFilter){
      //den
      if(bestTopCand.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	//num
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
	/*//Gen Matching
	if(!isData){
	  int match = std::min(partonMatch(bestTopCand, 0.6), constMatch(bestTopCand, 0.6));
	  if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	    if(match==-1 || match==0) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_0m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==1) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_1m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==2) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_2m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==3) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_3m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  }
	  if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	    if(match==-1 || match==0) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_0m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==1) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_1m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==2) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_2m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==3) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_3m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  }
	  if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	    if(match==-1 || match==0) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_0m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==1) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_1m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==2) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_2m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==3) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_3m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  }
	  if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	    if(match==-1 || match==0) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_0m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==1) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_1m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==2) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_2m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	    if(match==3) pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_3m,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	  }
	  }*/
      }
      //for events variables
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
      //num
      if(isTopPassDisc(Ntop,discThrLWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
      if(isTopPassDisc(Ntop,discThrMWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
      if(isTopPassDisc(Ntop,discThrTWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
      if(isTopPassDisc(Ntop,discThrXWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff_PDFDown);
	}
    }
    

    //Mistag Rate
    //QCD CS (data triggered with with ht)
    if(passNjet && passLepVeto && passQCDHT && passHtTrigger && Pass_EventFilter){
      //for best cand.                                                                                                                        
      //den                                                                                                                                   
      if(bestTopCand.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
        //num                                                                                                                                 
        if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
        }
        if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
        }
        if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
        }
	if(isbestCandTag && bestTopCand.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_BTUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_BTDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_PUUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_PUDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_PDFUp,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_PDFDown,bestTopCand.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
        }
      }
      //for events variables
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
      //num
      if(isTopPassDisc(Ntop,discThrLWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
	}
      if(isTopPassDisc(Ntop,discThrMWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
	}
      if(isTopPassDisc(Ntop,discThrTWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
	}
      if(isTopPassDisc(Ntop,discThrXWP))
	{
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_BTUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_BTDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_BTDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_PUUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_PUDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PUDown);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_PDFUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFUp);
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_PDFDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis_PDFDown);
	}

    }
    if(!isGenSys){ 
  // top tagger with JECUp
    //get output of tagger
    const TopTaggerResults* ttrJECUp = tr->getVar<const TopTaggerResults*>("ttrMVAJECUp");
    //Use result for top var
    vector<TopObject*> NtopJECUp = ttrJECUp->getTops();
    vector<TopObject> NtopCandJECUp = ttrJECUp->getTopCandidates();
    //std::cout<<"JECUp: "<<"Ntop: "<<NtopJECUp.size()<<" NtopCand: "<<NtopCandJECUp.size()<<std::endl;

    TopObject bestTopCandJECUp = GetBestTopCand(NtopCandJECUp);
    bool isbestCandTagJECUp = CandTag(bestTopCandJECUp, NtopJECUp);
    //eff
    if(passNjet && passMu && passBjet && passdPhi && passBLep && passMassBLep && passdPhiLep && passMtW && passHT && passLowMET && passExtraMuonVeto && passEleVeto && passMuTrigger && Pass_EventFilter){
      //den                                                                                                                                                             
      if(bestTopCandJECUp.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	//num                                                                                                                                                           
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      }
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
      //num                                                                                                                                                             
      if(isTopPassDisc(NtopJECUp,discThrLWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(NtopJECUp,discThrMWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(NtopJECUp,discThrTWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(NtopJECUp,discThrXWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
    }
    //mistagrate
    if(passNjet && passLepVeto && passQCDHT && passHtTrigger && Pass_EventFilter){
      //den                                                                                                                                                         
      if(bestTopCandJECUp.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	//num                                                                                                                                                           
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
	if(isbestCandTagJECUp && bestTopCandJECUp.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_JECUp,bestTopCandJECUp.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      }
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      //num                                                                                                                                                             
      if(isTopPassDisc(NtopJECUp,discThrLWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      if(isTopPassDisc(NtopJECUp,discThrMWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      if(isTopPassDisc(NtopJECUp,discThrTWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      if(isTopPassDisc(NtopJECUp,discThrXWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_JECUp,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}

    }

  // top tagger with JECDown
    //get output of tagger
    const TopTaggerResults* ttrJECDown = tr->getVar<const TopTaggerResults*>("ttrMVAJECDown");
    //Use result for top var
    vector<TopObject*> NtopJECDown = ttrJECDown->getTops();
    vector<TopObject> NtopCandJECDown = ttrJECDown->getTopCandidates();
    //std::cout<<"JECDown: "<<"Ntop: "<<NtopJECDown.size()<<" NtopCand: "<<NtopCandJECDown.size()<<std::endl;

    TopObject bestTopCandJECDown = GetBestTopCand(NtopCandJECDown);
    bool isbestCandTagJECDown = CandTag(bestTopCandJECDown, NtopJECDown);
    //eff
    if(passNjet && passMu && passBjet && passdPhi && passBLep && passMassBLep && passdPhiLep && passMtW && passHT && passLowMET && passExtraMuonVeto && passEleVeto && passMuTrigger && Pass_EventFilter){
      //den                                                                                                                                                             
      if(bestTopCandJECDown.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtEff_den_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	//num                                                                                                                                                           
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_LWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_MWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_TWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtEff_num_XWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      }
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetEff_den_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
      //num                                                                                                                                                             
      if(isTopPassDisc(NtopJECDown,discThrLWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_LWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(NtopJECDown,discThrMWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_MWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(NtopJECDown,discThrTWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_TWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
      if(isTopPassDisc(NtopJECDown,discThrXWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetEff_num_XWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_eff);
	}
    }
    //mistagrate
    if(passNjet && passLepVeto && passQCDHT && passHtTrigger && Pass_EventFilter){
      //den                                                                                                                                                         
      if(bestTopCandJECDown.p().M()>0.0){
	pUtility::FillDouble(myBaseHistgram.htopPtMis_den_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	//num                                                                                                                                                           
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrLWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_LWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrMWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_MWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrTWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_TWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
	if(isbestCandTagJECDown && bestTopCandJECDown.getDiscriminator()>=discThrXWP){
	  pUtility::FillDouble(myBaseHistgram.htopPtMis_num_XWP_JECDown,bestTopCandJECDown.p().Pt(),Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      }
      //den
      pUtility::FillDouble(myBaseHistgram.hNjetMis_den_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
      //num                                                                                                                                                             
      if(isTopPassDisc(NtopJECDown,discThrLWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_LWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      if(isTopPassDisc(NtopJECDown,discThrMWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_MWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      if(isTopPassDisc(NtopJECDown,discThrTWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_TWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}
      if(isTopPassDisc(NtopJECDown,discThrXWP))
        {
	  pUtility::FillDouble(myBaseHistgram.hNjetMis_num_XWP_JECDown,cntNJetsPt30Eta24,Lumiscale * EventWeight*tot_MCWgt_mis);
	}

    }
   }
  }//event loop
  (myBaseHistgram.oFile)->Write();
  return 0;
}
