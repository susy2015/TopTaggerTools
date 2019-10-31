#include <iostream>
#include <string>
#include <vector>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TAxis.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TProfile.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"

#include "ResTagger.h"

using namespace std;

double Lumiscale = 1.0;
double EventWeight = 1.0;

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;

  TH1D *hrecotopPt_LWP;
  TH1D *hrecotopPt_MWP;
  TH1D *hrecotopPt_TWP;
  TH1D *hrecotopPt_XWP;
  TH1D *hrecotopPt;

  TH1D *hgentopPt;
  TH1D *hstopMass;
  TH1D *hneutralinoMass;

  TH1D *htopPtmis_den;
  TH1D *htopPtmisB_den;
  TH1D *htopPtmisNoB_den;
  TH1D *htopPtmis_num_LWP;
  TH1D *htopPtmisB_num_LWP;
  TH1D *htopPtmisNoB_num_LWP;
  TH1D *htopPtmis_num_MWP;
  TH1D *htopPtmisB_num_MWP;
  TH1D *htopPtmisNoB_num_MWP;
  TH1D *htopPtmis_num_TWP;
  TH1D *htopPtmisB_num_TWP;
  TH1D *htopPtmisNoB_num_TWP;
  TH1D *htopPtmis_num_XWP;
  TH1D *htopPtmisB_num_XWP;
  TH1D *htopPtmisNoB_num_XWP;

  TH1D *hNjetmis_den;
  TH1D *hNjetmis_num_LWP;
  TH1D *hNjetmis_num_MWP;
  TH1D *hNjetmis_num_TWP;
  TH1D *hNjetmis_num_XWP;

};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_Eff"+index+".root";
  oFile = new TFile(filename, "recreate");

  hrecotopPt_LWP = new TH1D("hrecotopPt_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_LWP->Sumw2();
  hrecotopPt_MWP = new TH1D("hrecotopPt_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_MWP->Sumw2();
  hrecotopPt_TWP = new TH1D("hrecotopPt_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_TWP->Sumw2();
  hrecotopPt_XWP = new TH1D("hrecotopPt_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_XWP->Sumw2();
  hrecotopPt = new TH1D("hrecotopPt", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt->Sumw2();
  hgentopPt = new TH1D("hgentopPt", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hgentopPt->Sumw2();
  hstopMass = new TH1D("hstopMass", "Stop Mass;M[GeV];Event", 40, 500, 2500);
  hstopMass->Sumw2();
  hneutralinoMass = new TH1D("hneutralinoMass", "Neutralino Mass;M[GeV];Event", 24, 50, 1250);
  hneutralinoMass->Sumw2();

  htopPtmis_den = new TH1D("htopPtmis_den", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmis_den->Sumw2();
  htopPtmis_num_LWP = new TH1D("htopPtmis_num_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmis_num_LWP->Sumw2();
  htopPtmis_num_MWP = new TH1D("htopPtmis_num_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmis_num_MWP->Sumw2();
  htopPtmis_num_TWP = new TH1D("htopPtmis_num_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmis_num_TWP->Sumw2();
  htopPtmis_num_XWP = new TH1D("htopPtmis_num_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmis_num_XWP->Sumw2();

  htopPtmisB_den = new TH1D("htopPtmisB_den", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisB_den->Sumw2();
  htopPtmisB_num_LWP = new TH1D("htopPtmisB_num_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisB_num_LWP->Sumw2();
  htopPtmisB_num_MWP = new TH1D("htopPtmisB_num_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisB_num_MWP->Sumw2();
  htopPtmisB_num_TWP = new TH1D("htopPtmisB_num_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisB_num_TWP->Sumw2();
  htopPtmisB_num_XWP = new TH1D("htopPtmisB_num_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisB_num_XWP->Sumw2();

  htopPtmisNoB_den = new TH1D("htopPtmisNoB_den", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisNoB_den->Sumw2();
  htopPtmisNoB_num_LWP = new TH1D("htopPtmisNoB_num_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisNoB_num_LWP->Sumw2();
  htopPtmisNoB_num_MWP = new TH1D("htopPtmisNoB_num_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisNoB_num_MWP->Sumw2();
  htopPtmisNoB_num_TWP = new TH1D("htopPtmisNoB_num_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisNoB_num_TWP->Sumw2();
  htopPtmisNoB_num_XWP = new TH1D("htopPtmisNoB_num_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtmisNoB_num_XWP->Sumw2();

  hNjetmis_den = new TH1D("hNjetmis_den", "N_{jets};N_{jets};Event",16, 4, 20);
  hNjetmis_den->Sumw2();
  hNjetmis_num_LWP = new TH1D("hNjetmis_num_LWP", "N_{jets};N_{jets};Event",16, 4, 20);
  hNjetmis_num_LWP->Sumw2();
  hNjetmis_num_MWP = new TH1D("hNjetmis_num_MWP", "N_{jets};N_{jets};Event",16, 4, 20);
  hNjetmis_num_MWP->Sumw2();
  hNjetmis_num_TWP = new TH1D("hNjetmis_num_TWP", "N_{jets};N_{jets};Event",16, 4, 20);
  hNjetmis_num_TWP->Sumw2();
  hNjetmis_num_XWP = new TH1D("hNjetmis_num_XWP", "N_{jets};N_{jets};Event",16, 4, 20);
  hNjetmis_num_XWP->Sumw2();
  
}

bool FillChain(TChain* &chain, const char *subsample, const bool iscondor, const int& startfile, const int& filerun){
  AnaSamples::SampleSet        allSamples = AnaSamples::SampleSet("sampleSets_PostProcessed_2016.cfg", iscondor, AnaSamples::luminosity_muon_2016);
  AnaSamples::SampleCollection allCollections("sampleCollections_2016.cfg", allSamples);
  bool find = false;  
  TString subsamplename(subsample);
  chain = new TChain(allSamples[subsample].treePath.c_str());
  if(allSamples[subsample] != allSamples.null())
    {
      allSamples[subsample].addFilesToChain(chain, startfile, filerun);
      find = true;
      Lumiscale = allSamples[subsample].getWeight();  
    }
  return find;
}

std::vector<int> GetMuIdx(const std::vector<TLorentzVector> &muonsLVec, const std::vector<float> &muonsMiniIso, const std::vector<float> &muonsMtw, const std::vector<unsigned char> &muonsFlagIDVec, const AnaConsts::IsoAccRec& muonsMiniIsoArr){
  std::vector<int> Idx;
  for(int im=0; im<muonsLVec.size(); im++){
    if(AnaFunctions::passMuon(muonsLVec[im], muonsMiniIso[im], muonsMtw[im], muonsFlagIDVec[im], muonsMiniIsoArr))
      {
        Idx.push_back(im);
      }
  }
  return Idx;

}

TLorentzVector GetMuLVec(const std::vector<TLorentzVector> &muonsLVec, const std::vector<int> &muonsIdx, int &id, const float &ptCut){
  TLorentzVector passLep(0.0, 0.0, 0.0, 0.0);
  for(int im=0; im<muonsIdx.size(); im++){
    if(muonsLVec[muonsIdx[im]].Pt()>ptCut) {
      passLep = muonsLVec[muonsIdx[im]];
      id = muonsIdx[im];
      break;
    }
  }
  return passLep;
}


std::vector<TLorentzVector> GetLep(const std::vector<TLorentzVector> &muonsLVec, const std::vector<int> &muonsIdx, const int &id){
  std::vector<TLorentzVector> passLep;
  for(int im=0; im<muonsIdx.size(); im++){
    if(muonsIdx[im]==id) continue;
    passLep.push_back(muonsLVec[muonsIdx[im]]);
  }
  return passLep;
}

bool CheckDiLep(const std::vector<TLorentzVector> &muonsLVec, const std::vector<float> &muonsMiniIso, const std::vector<float> &muonsMtw, const std::vector<unsigned char> &muonsFlagIDVec, const AnaConsts::IsoAccRec& muonsMiniIsoArr, const std::vector<TLorentzVector> &elesLVec, const std::vector<float> &elesMiniIso, const std::vector<float> &elesMtw, const std::vector<unsigned char> &elesFlagIDVec, const AnaConsts::IsoAccRec& elesMiniIsoArr, const std::vector<int> &muonsChrg, const std::vector<int> &elesChrg, const float &ptCut1, const float &ptCut2){
  bool check = false;
  for(int im=0; im<muonsLVec.size(); im++){
    if(AnaFunctions::passMuon(muonsLVec[im], muonsMiniIso[im], muonsMtw[im], muonsFlagIDVec[im], muonsMiniIsoArr) && (muonsLVec[im].Pt()>ptCut1))
      {
	for(int ie=0; ie<elesLVec.size(); ie++){
	  if(AnaFunctions::passElectron(elesLVec[ie], elesMiniIso[ie], elesMtw[ie], elesFlagIDVec[ie], elesMiniIsoArr) && (elesLVec[ie].Pt()>ptCut2))
	    {
	      if((muonsChrg[im]*elesChrg[ie])<0){check=true; break;}
	    }
	}
      }
  }
  return check;
}

TopObject GetBestTopCand(std::vector<TopObject> &NtopCand){
  TopObject bestCand;
  double bestTopMass = -9999.9;
  for(auto cand:NtopCand){
    if(fabs(cand.p().M() - 173.5) < fabs(bestTopMass - 173.5) && cand.getNConstituents() == 3)
      {
	bestTopMass = cand.p().M();
	bestCand = cand;
      }
  }
  return bestCand;
}

/*TopObject GetBestTopCand(std::vector<TopObject> &NtopCand){
  TopObject bestCand;
  double bestTopDisc = -9999.9;
  for(auto cand:NtopCand){
    if(cand.getDiscriminator() > bestTopDisc)
      {
	bestTopDisc = cand.getDiscriminator();
	bestCand = cand;
      }
  }
  return bestCand;
  }*/

bool CandTag(TopObject bestCand, std::vector<TopObject*> tops){
  bool bestCandTag = false;
  TopObject* bestCandptr = &bestCand;
  for(const auto& top : tops)
    {
      if(top->p() == bestCandptr->p())
	{
	  bestCandTag = true;
	  break;
	}
    }
  return bestCandTag;
}

bool isTopPassDisc(std::vector<TopObject*> tops, const double disc){
  bool ispass = false;
  for(const auto& top : tops)
    {
      if(top->getDiscriminator() >= disc)
        {
          ispass = true;
          break;
        }
    }

  return ispass;
}

int partonMatch(TopObject bestCand, double dRMax){
  int match = -1;
  double bestMatchDR = 999.9;
  const TLorentzVector* bestMatch = nullptr;
  for(const auto& genTop : bestCand.getGenTopMatches()){
    double deltaR = ROOT::Math::VectorUtil::DeltaR(bestCand.p(), *genTop.first);
    if(deltaR < bestMatchDR)
      {
	bestMatchDR = deltaR;
	bestMatch = genTop.first;
      }
  }
  //if(bestMatchDR < dRMax) match = bestCand.getGenTopMatches()[bestMatch].size();
  if(bestMatch != nullptr) match = bestCand.getGenTopMatches().at(bestMatch).size();
  return match;
}

int constMatch(TopObject bestCand, double dRMax){
  int match = -1;
  double bestMatchDR = 999.9;
  const TLorentzVector* bestMatch = nullptr;
  int nConstMatch = -1;                                                                                                                    
  for(const auto& genTop : bestCand.getGenTopMatches()){
    nConstMatch = 0;
    for(const auto& constituent : bestCand.getConstituents()){                                                                                                    
      const auto& genTopIter = constituent->getGenMatches().find(genTop.first);                                                                                      
      if(genTopIter != constituent->getGenMatches().end()) ++nConstMatch; 
    }
    double deltaR = ROOT::Math::VectorUtil::DeltaR(bestCand.p(), *genTop.first);
    if(deltaR < bestMatchDR)
      {
	bestMatchDR = deltaR;
	bestMatch = genTop.first;
      }
  }
  //if(bestMatchDR < dRMax) match = bestCand.getGenTopMatches()[bestMatch].size();
  if(bestMatch != nullptr) match = nConstMatch;
  return match;
}

bool PassHEM(const std::vector<TLorentzVector> &jetsLVec, const float &etalow, const float &etahigh, const float &philow, const float &phihigh, const float &pt){
  bool hem = true;
  for(int i=0; i<jetsLVec.size(); i++){
    if((jetsLVec[i].Eta() >= etalow && jetsLVec[i].Eta() <= etahigh) && (jetsLVec[i].Phi() >= philow && jetsLVec[i].Phi() <= phihigh) && (jetsLVec[i].Pt() > pt)) 
      {
	hem = false;
	break;
      }

  }
  return hem;
}

bool CheckBconst(TopObject bestCand, const float &cut){
  bool check = false;
  for(const auto& constituent : bestCand.getConstituents()){
    if(constituent->getBTagDisc()>cut && constituent->p().Pt()>20){
      check = true;
      break;
    }
  }
  return check;
}
