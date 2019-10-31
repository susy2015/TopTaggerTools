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
#include "TagVar.h"

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 4 arguments "<<"SubsampleName"<<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Tagvar TTbarInc 1000 0 1" << std::endl;
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
  //myBaseHistgram.BookHistgram(subsamplename, startfile);
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


    const vector<TLorentzVector> &jetsLVec = tr->getVec_LVFromNano<float>("Jet");
    
    float met=tr->getVar<float>("MET_pt");
    float metphi=tr->getVar<float>("MET_phi");
    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0.0, metphi, 0.0);


    float EvtWt = isData? 1 : tr->getVar<float>("Generator_weight");
    //std::cout<<"EvtWt"<<EvtWt<<std::endl;

    const unsigned int run = tr->getVar<unsigned int>("run");

    if(!isData){
      const vector<TLorentzVector> &genDecayLVec = tr->getVec_LVFromNano<float>("GenPart");
      const vector<int> &genDecayPdgIdVec = tr->getVec<int>("GenPart_pdgId");
      const vector<int> &genDecayMomIdxVec = tr->getVec<int>("GenPart_genPartIdxMother");
    }


    int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);

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
    //const bool Pass_EventFilter = tr->getVar<bool>("Pass_EventFilter");

    //MC weights
    float BTagWgt = 1.0;
    float PileUpWgt = 1.0;
    float TTbarWgt = 1.0;

    float tot_MCWgt_eff = 1.0;
    float tot_MCWgt_mis = 1.0;

    if(!isData){
      //pileUp weight
      //PileUpWgt = tr->getVar<float>("puWeight");
      //BTag weight
      //BTagWgt = tr->getVar<float>("BTagWeight");
      //TTbar weight
      bool doTTbarwgt = sample.Contains("TTbar")? true:false;
      if(doTTbarwgt){
	//TTbarWgt = tr->getVar<float>("ISRWeight");
      }
      tot_MCWgt_eff = PileUpWgt * BTagWgt * TTbarWgt;
      tot_MCWgt_mis = PileUpWgt * BTagWgt * TTbarWgt;
    }
    
    EventWeight = EvtWt>0? 1:-1;

    //cut
    bool passNjet(false);
    if(cntNJetsPt30Eta24>=4) passNjet = true;
    if(!passNjet) continue;//to make sure resolved tagger has at least 3 jets    
    // top tagger
    //get output of tagger
    const TopTaggerResults* ttr = tr->getVar<const TopTaggerResults*>("ttrMVA");
    //Use result for top var
    vector<TopObject*> Ntop = ttr->getTops();
    vector<TopObject> NtopCand = ttr->getTopCandidates();
    //cout<<"Ntop:  "<<Ntop.size()<<endl;
    //if(Ntop.size()==2)cout<<"Ntop[0]->getDiscriminator(): "<<Ntop[0]->getDiscriminator()<<" Ntop[1]->getDiscriminator(): "<<Ntop[1]->getDiscriminator()<<endl;
    if(Ntop.size()){
      map<string, float> Inputvar = Ntop[0]->getMVAInputs();
      //cout<<"TotalInputs:  "<<Inputvar.size()<<endl;
      for(auto var:Inputvar){
	//cout<<var.first<<" : "<<var.second<<endl;	
	TString name = var.first;
	float val = var.second;
	float nval = 1000;
	if(name.Contains("DeepCSV") || name.Contains("EnergyFraction") || name.Contains("qgAxis")) nval = 1;
	if(name.Contains("qgPtD") || name.Contains("Theta")) nval = 10;
	if(name.Contains("qgMult") || name.Contains("Multiplicity")) nval = 100;
	if(VarMap.size()<Inputvar.size()){
	  VarMap[name] = new TH1D(name, name+";InputVar;Events",50,val-val,nval);
	  VarMap[name]->Sumw2();
	  VarMap[name]->Fill(val, Lumiscale * EventWeight);
	}
	else  VarMap[name]->Fill(val, Lumiscale * EventWeight);
      }
    }
  }//event loop
  // for(auto var:VarMap){
  //myBaseHistgram.hvars.push_back(static_cast<TH1*>(var.second));
  //}
  //cout<<"VarMap.size(): "<<VarMap.size()<<endl;
  myBaseHistgram.BookHistgram(subsamplename, startfile, VarMap);
  (myBaseHistgram.oFile)->Write();
  return 0;
}
