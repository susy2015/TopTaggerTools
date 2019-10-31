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
map<TString, TH1*> VarMap;

class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&, const map<TString, TH1*>&);
  TFile *oFile;
  vector<TH1*> hvars;

};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun, const map<TString, TH1*>& varmap)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_TagVar"+index+".root";
  oFile = new TFile(filename, "recreate");
  for(auto var:varmap){
    hvars.push_back(static_cast<TH1D*>((var.second)->Clone(var.first)));
  }
}

bool FillChain(TChain* &chain, const char *subsample, const bool iscondor, const int& startfile, const int& filerun){
  AnaSamples::SampleSet        allSamples = AnaSamples::SampleSet("sampleSets_PostProcessed_2018.cfg", iscondor, AnaSamples::luminosity_muon_2018);
  AnaSamples::SampleCollection allCollections("sampleCollections_2018.cfg", allSamples);
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


