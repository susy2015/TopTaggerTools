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
  TH1D *hNtop;
  TH1D *hNtopCand;
  TH1D *htopCandDisc;
  TH1D *hNjet_Eff;
  TH1D *hmet_Eff;
  TH1D *hht_Eff;
  TH1D *hBestCandDisc_Eff;
  TH1D *hNjet_Mis;
  TH1D *hmet_Mis;
  TH1D *hht_Mis;
  TH1D *hBestCandDisc_Mis;
  TH1D *hNjet_2l;
  TH1D *hmet_2l;
  TH1D *hht_2l;
  TH1D *hBestCandDisc_2l;

  TH1D *hcutflow_eff;
  TH1D *hcutflow_mis;

  TH1D *htopPtEff_den;
  TH1D *htopPtEff_num_LWP;
  TH1D *htopPtEff_num_MWP;
  TH1D *htopPtEff_num_TWP;
  TH1D *htopPtEff_num_XWP;
  TH1D *hNjetEff_den;
  TH1D *hNjetEff_num_LWP;
  TH1D *hNjetEff_num_MWP;
  TH1D *hNjetEff_num_TWP;
  TH1D *hNjetEff_num_XWP;
  TH1D *hmetEff_den;
  TH1D *hmetEff_num_LWP;
  TH1D *hmetEff_num_MWP;
  TH1D *hmetEff_num_TWP;
  TH1D *hmetEff_num_XWP;
  TH1D *hhtEff_den;
  TH1D *hhtEff_num_LWP;
  TH1D *hhtEff_num_MWP;
  TH1D *hhtEff_num_TWP;
  TH1D *hhtEff_num_XWP;
  TH2D *htopPtNjetEff_den;
  TH2D *htopPtNjetEff_num_LWP;
  TH2D *htopPtNjetEff_num_MWP;
  TH2D *htopPtNjetEff_num_TWP;
  TH2D *htopPtNjetEff_num_XWP;

  TH1D *htopPtMis_den;
  TH1D *htopPtMis_num_LWP;
  TH1D *htopPtMis_num_MWP;
  TH1D *htopPtMis_num_TWP;
  TH1D *htopPtMis_num_XWP;
  TH1D *hNjetMis_den;
  TH1D *hNjetMis_num_LWP;
  TH1D *hNjetMis_num_MWP;
  TH1D *hNjetMis_num_TWP;
  TH1D *hNjetMis_num_XWP;
  TH1D *hmetMis_den;
  TH1D *hmetMis_num_LWP;
  TH1D *hmetMis_num_MWP;
  TH1D *hmetMis_num_TWP;
  TH1D *hmetMis_num_XWP;
  TH1D *hhtMis_den;
  TH1D *hhtMis_num_LWP;
  TH1D *hhtMis_num_MWP;
  TH1D *hhtMis_num_TWP;
  TH1D *hhtMis_num_XWP;
  TH2D *htopPtNjetMis_den;
  TH2D *htopPtNjetMis_num_LWP;
  TH2D *htopPtNjetMis_num_MWP;
  TH2D *htopPtNjetMis_num_TWP;
  TH2D *htopPtNjetMis_num_XWP;

  TH1D *htopPt2l_den;
  TH1D *htopPt2l_num_LWP;
  TH1D *htopPt2l_num_MWP;
  TH1D *htopPt2l_num_TWP;
  TH1D *htopPt2l_num_XWP;
  TH1D *hNjet2l_den;
  TH1D *hNjet2l_num_LWP;
  TH1D *hNjet2l_num_MWP;
  TH1D *hNjet2l_num_TWP;
  TH1D *hNjet2l_num_XWP;  

  TH2D *hNjetDisc_Eff;
  TH2D *hNjetDisc_Mis;

  TH1D *hrecotopPt_LWP;
  TH1D *hrecotopPt_0m_LWP;
  TH1D *hrecotopPt_1m_LWP;
  TH1D *hrecotopPt_2m_LWP;
  TH1D *hrecotopPt_3m_LWP;
  TH1D *hrecotopPt_MWP;
  TH1D *hrecotopPt_0m_MWP;
  TH1D *hrecotopPt_1m_MWP;
  TH1D *hrecotopPt_2m_MWP;
  TH1D *hrecotopPt_3m_MWP;
  TH1D *hrecotopPt_TWP;
  TH1D *hrecotopPt_0m_TWP;
  TH1D *hrecotopPt_1m_TWP;
  TH1D *hrecotopPt_2m_TWP;
  TH1D *hrecotopPt_3m_TWP;
  TH1D *hrecotopPt_XWP;
  TH1D *hrecotopPt_0m_XWP;
  TH1D *hrecotopPt_1m_XWP;
  TH1D *hrecotopPt_2m_XWP;
  TH1D *hrecotopPt_3m_XWP;
  TH1D *hrecotopPt;
  TH1D *hrecotopPt_0m;
  TH1D *hrecotopPt_1m;
  TH1D *hrecotopPt_2m;
  TH1D *hrecotopPt_3m;

  TH1D *hrecotopPt_Mis_LWP;
  TH1D *hrecotopPt_Mis_0m_LWP;
  TH1D *hrecotopPt_Mis_1m_LWP;
  TH1D *hrecotopPt_Mis_2m_LWP;
  TH1D *hrecotopPt_Mis_3m_LWP;
  TH1D *hrecotopPt_Mis_MWP;
  TH1D *hrecotopPt_Mis_0m_MWP;
  TH1D *hrecotopPt_Mis_1m_MWP;
  TH1D *hrecotopPt_Mis_2m_MWP;
  TH1D *hrecotopPt_Mis_3m_MWP;
  TH1D *hrecotopPt_Mis_TWP;
  TH1D *hrecotopPt_Mis_0m_TWP;
  TH1D *hrecotopPt_Mis_1m_TWP;
  TH1D *hrecotopPt_Mis_2m_TWP;
  TH1D *hrecotopPt_Mis_3m_TWP;
  TH1D *hrecotopPt_Mis_XWP;
  TH1D *hrecotopPt_Mis_0m_XWP;
  TH1D *hrecotopPt_Mis_1m_XWP;
  TH1D *hrecotopPt_Mis_2m_XWP;
  TH1D *hrecotopPt_Mis_3m_XWP;
  TH1D *hrecotopPt_Mis;
  TH1D *hrecotopPt_Mis_0m;
  TH1D *hrecotopPt_Mis_1m;
  TH1D *hrecotopPt_Mis_2m;
  TH1D *hrecotopPt_Mis_3m;

  TH1D *hNjetEff_LWP;
  TH1D *hNjetEff_0m_LWP;
  TH1D *hNjetEff_1m_LWP;
  TH1D *hNjetEff_2m_LWP;
  TH1D *hNjetEff_3m_LWP;
  TH1D *hNjetEff_MWP;
  TH1D *hNjetEff_0m_MWP;
  TH1D *hNjetEff_1m_MWP;
  TH1D *hNjetEff_2m_MWP;
  TH1D *hNjetEff_3m_MWP;
  TH1D *hNjetEff_TWP;
  TH1D *hNjetEff_0m_TWP;
  TH1D *hNjetEff_1m_TWP;
  TH1D *hNjetEff_2m_TWP;
  TH1D *hNjetEff_3m_TWP;
  TH1D *hNjetEff_XWP;
  TH1D *hNjetEff_0m_XWP;
  TH1D *hNjetEff_1m_XWP;
  TH1D *hNjetEff_2m_XWP;
  TH1D *hNjetEff_3m_XWP;
  TH1D *hNjetEff;
  TH1D *hNjetEff_0m;
  TH1D *hNjetEff_1m;
  TH1D *hNjetEff_2m;
  TH1D *hNjetEff_3m;

  //Rebin                                                                                                              
  double ptbins[6] = {0.0, 250.0, 300.0, 350.0, 450.0, 1000.0};
  double jetbins[6]= {4.0, 5.0, 6.0, 7.0, 8.0, 15.0};
  int ptarr = sizeof(ptbins)/sizeof(ptbins[0]);
  int ptbin = ptarr-1;
  int jetarr = sizeof(jetbins)/sizeof(jetbins[0]);
  int jetbin = jetarr-1;
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_DataMC"+index+".root";
  oFile = new TFile(filename, "recreate");
 
  hNtop = new TH1D("hNtop", "No. of top;N_{top};Event", 10, 0, 10);
  hNtop->Sumw2();
  hNtopCand = new TH1D("hNtopCand", "No. of top candidate;N_{topCand};Event", 10, 0, 10);
  hNtopCand->Sumw2();
  htopCandDisc = new TH1D("htopCandDisc", "Cand. discriminator;Disc.;Event", 100, 0, 1);
  htopCandDisc->Sumw2();
  hNjet_Eff = new TH1D("hNjet_Eff","N_{jets};N_{jet};Event",15,0,15);
  hNjet_Eff->Sumw2();
  hmet_Eff = new TH1D("hmet_Eff","met;p_{T}^{miss}[GeV];Event",20, 0, 1000);
  hmet_Eff->Sumw2();
  hht_Eff = new TH1D("hht_Eff","ht;H_{T}[GeV];Event",60, 0, 3000);
  hht_Eff->Sumw2();
  hBestCandDisc_Eff = new TH1D("hBestCandDisc_Eff", "Cand. discriminator;Disc.;Event", 100, 0, 1);
  hBestCandDisc_Eff->Sumw2();
  hNjet_Mis = new TH1D("hNjet_Mis","N_{jets};N_{jet};Event",15,0,15);
  hNjet_Mis->Sumw2();
  hmet_Mis = new TH1D("hmet_Mis","met;p_{T}^{miss}[GeV];Event",20, 0, 1000);
  hmet_Mis->Sumw2();
  hht_Mis = new TH1D("hht_Mis","ht;H_{T}[GeV];Event",60, 0, 3000);
  hht_Mis->Sumw2();
  hBestCandDisc_Mis = new TH1D("hBestCandDisc_Mis", "Cand. discriminator;Disc.;Event", 100, 0, 1);
  hBestCandDisc_Mis->Sumw2();
  hNjet_2l = new TH1D("hNjet_2l","N_{jets};N_{jet};Event",15,0,15);
  hNjet_2l->Sumw2();
  hmet_2l = new TH1D("hmet_2l","met;p_{T}^{miss}[GeV];Event",20, 0, 1000);
  hmet_2l->Sumw2();
  hht_2l = new TH1D("hht_2l","ht;H_{T}[GeV];Event",60, 0, 3000);
  hht_2l->Sumw2();
  hBestCandDisc_2l = new TH1D("hBestCandDisc_2l", "Cand. discriminator;Disc.;Event", 100, 0, 1);
  hBestCandDisc_2l->Sumw2();

  hcutflow_eff = new TH1D("hcutflow_eff", "cut flow table", 20, 0, 20);
  hcutflow_eff->SetCanExtend(TH1::kXaxis);                                                       
  hcutflow_eff->Sumw2();
  hcutflow_mis = new TH1D("hcutflow_mis", "cut flow table", 20, 0, 20);
  hcutflow_mis->SetCanExtend(TH1::kXaxis);                                           
  hcutflow_mis->Sumw2();


  htopPtEff_den = new TH1D("htopPtEff_den", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den->Sumw2();
  htopPtEff_num_LWP = new TH1D("htopPtEff_num_LWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP->Sumw2();
  htopPtEff_num_MWP = new TH1D("htopPtEff_num_MWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP->Sumw2();
  htopPtEff_num_TWP = new TH1D("htopPtEff_num_TWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP->Sumw2();
  htopPtEff_num_XWP = new TH1D("htopPtEff_num_XWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP->Sumw2();
  hNjetEff_den = new TH1D("hNjetEff_den","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den->Sumw2();
  hNjetEff_num_LWP = new TH1D("hNjetEff_num_LWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP->Sumw2();
  hNjetEff_num_MWP = new TH1D("hNjetEff_num_MWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP->Sumw2();
  hNjetEff_num_TWP = new TH1D("hNjetEff_num_TWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP->Sumw2();
  hNjetEff_num_XWP = new TH1D("hNjetEff_num_XWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP->Sumw2();
  hmetEff_den = new TH1D("hmetEff_den","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetEff_den->Sumw2();
  hmetEff_num_LWP = new TH1D("hmetEff_num_LWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetEff_num_LWP->Sumw2();
  hmetEff_num_MWP = new TH1D("hmetEff_num_MWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetEff_num_MWP->Sumw2();
  hmetEff_num_TWP = new TH1D("hmetEff_num_TWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetEff_num_TWP->Sumw2();
  hmetEff_num_XWP = new TH1D("hmetEff_num_XWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetEff_num_XWP->Sumw2();
  hhtEff_den = new TH1D("hhtEff_den","ht;H_{T}[GeV];Event",60,0,3000);
  hhtEff_den->Sumw2();
  hhtEff_num_LWP = new TH1D("hhtEff_num_LWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtEff_num_LWP->Sumw2();
  hhtEff_num_MWP = new TH1D("hhtEff_num_MWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtEff_num_MWP->Sumw2();
  hhtEff_num_TWP = new TH1D("hhtEff_num_TWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtEff_num_TWP->Sumw2();
  hhtEff_num_XWP = new TH1D("hhtEff_num_XWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtEff_num_XWP->Sumw2();
  htopPtNjetEff_den = new TH2D("htopPtNjetEff_den", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetEff_den->Sumw2();
  htopPtNjetEff_num_LWP = new TH2D("htopPtNjetEff_num_LWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetEff_num_LWP->Sumw2();
  htopPtNjetEff_num_MWP = new TH2D("htopPtNjetEff_num_MWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetEff_num_MWP->Sumw2();
  htopPtNjetEff_num_TWP = new TH2D("htopPtNjetEff_num_TWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetEff_num_TWP->Sumw2();
  htopPtNjetEff_num_XWP = new TH2D("htopPtNjetEff_num_XWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetEff_num_XWP->Sumw2();

  htopPtMis_den = new TH1D("htopPtMis_den", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den->Sumw2();
  htopPtMis_num_LWP = new TH1D("htopPtMis_num_LWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP->Sumw2();
  htopPtMis_num_MWP = new TH1D("htopPtMis_num_MWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP->Sumw2();
  htopPtMis_num_TWP = new TH1D("htopPtMis_num_TWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP->Sumw2();
  htopPtMis_num_XWP = new TH1D("htopPtMis_num_XWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP->Sumw2();
  hNjetMis_den = new TH1D("hNjetMis_den","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den->Sumw2();
  hNjetMis_num_LWP = new TH1D("hNjetMis_num_LWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP->Sumw2();
  hNjetMis_num_MWP = new TH1D("hNjetMis_num_MWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP->Sumw2();
  hNjetMis_num_TWP = new TH1D("hNjetMis_num_TWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP->Sumw2();
  hNjetMis_num_XWP = new TH1D("hNjetMis_num_XWP","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP->Sumw2();
  hmetMis_den = new TH1D("hmetMis_den","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetMis_den->Sumw2();
  hmetMis_num_LWP = new TH1D("hmetMis_num_LWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetMis_num_LWP->Sumw2();
  hmetMis_num_MWP = new TH1D("hmetMis_num_MWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetMis_num_MWP->Sumw2();
  hmetMis_num_TWP = new TH1D("hmetMis_num_TWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetMis_num_TWP->Sumw2();
  hmetMis_num_XWP = new TH1D("hmetMis_num_XWP","met;p_{T}^{miss}[GeV];Event",20,0,1000);
  hmetMis_num_XWP->Sumw2();
  hhtMis_den = new TH1D("hhtMis_den","ht;H_{T}[GeV];Event",60,0,3000);
  hhtMis_den->Sumw2();
  hhtMis_num_LWP = new TH1D("hhtMis_num_LWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtMis_num_LWP->Sumw2();
  hhtMis_num_MWP = new TH1D("hhtMis_num_MWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtMis_num_MWP->Sumw2();
  hhtMis_num_TWP = new TH1D("hhtMis_num_TWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtMis_num_TWP->Sumw2();
  hhtMis_num_XWP = new TH1D("hhtMis_num_XWP","ht;H_{T}[GeV];Event",60,0,3000);
  hhtMis_num_XWP->Sumw2();
  htopPtNjetMis_den = new TH2D("htopPtNjetMis_den", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetMis_den->Sumw2();
  htopPtNjetMis_num_LWP = new TH2D("htopPtNjetMis_num_LWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetMis_num_LWP->Sumw2();
  htopPtNjetMis_num_MWP = new TH2D("htopPtNjetMis_num_MWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetMis_num_MWP->Sumw2();
  htopPtNjetMis_num_TWP = new TH2D("htopPtNjetMis_num_TWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetMis_num_TWP->Sumw2();
  htopPtNjetMis_num_XWP = new TH2D("htopPtNjetMis_num_XWP", "Efficiency;p_{T}[GeV];N_{jet}",ptbin,ptbins,jetbin,jetbins);
  htopPtNjetMis_num_XWP->Sumw2();

  htopPt2l_den = new TH1D("htopPt2l_den", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPt2l_den->Sumw2();
  htopPt2l_num_LWP = new TH1D("htopPt2l_num_LWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPt2l_num_LWP->Sumw2();
  htopPt2l_num_MWP = new TH1D("htopPt2l_num_MWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPt2l_num_MWP->Sumw2();
  htopPt2l_num_TWP = new TH1D("htopPt2l_num_TWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPt2l_num_TWP->Sumw2();
  htopPt2l_num_XWP = new TH1D("htopPt2l_num_XWP", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPt2l_num_XWP->Sumw2();
  hNjet2l_den = new TH1D("hNjet2l_den","N_{jets};N_{jet};Event",15,0,15);
  hNjet2l_den->Sumw2();
  hNjet2l_num_LWP = new TH1D("hNjet2l_num_LWP","N_{jets};N_{jet};Event",15,0,15);
  hNjet2l_num_LWP->Sumw2();
  hNjet2l_num_MWP = new TH1D("hNjet2l_num_MWP","N_{jets};N_{jet};Event",15,0,15);
  hNjet2l_num_MWP->Sumw2();
  hNjet2l_num_TWP = new TH1D("hNjet2l_num_TWP","N_{jets};N_{jet};Event",15,0,15);
  hNjet2l_num_TWP->Sumw2();
  hNjet2l_num_XWP = new TH1D("hNjet2l_num_XWP","N_{jets};N_{jet};Event",15,0,15);
  hNjet2l_num_XWP->Sumw2();

  hNjetDisc_Eff = new TH2D("hNjetDisc_Eff","N_{jets} vs Disc.;N_{jet};Disc.",10,4,14, 10, 0.5, 1.0);
  hNjetDisc_Eff->Sumw2();
  hNjetDisc_Mis = new TH2D("hNjetDisc_Mis","N_{jets} vs Disc.;N_{jet};Disc.",10,4,14, 10, 0.5, 1.0);
  hNjetDisc_Mis->Sumw2();

  hrecotopPt_LWP = new TH1D("hrecotopPt_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_LWP->Sumw2();
  hrecotopPt_0m_LWP = new TH1D("hrecotopPt_0m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_0m_LWP->Sumw2();
  hrecotopPt_1m_LWP = new TH1D("hrecotopPt_1m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_1m_LWP->Sumw2();
  hrecotopPt_2m_LWP = new TH1D("hrecotopPt_2m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_2m_LWP->Sumw2();
  hrecotopPt_3m_LWP = new TH1D("hrecotopPt_3m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_3m_LWP->Sumw2();
  hrecotopPt_MWP = new TH1D("hrecotopPt_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_MWP->Sumw2();
  hrecotopPt_0m_MWP = new TH1D("hrecotopPt_0m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_0m_MWP->Sumw2();
  hrecotopPt_1m_MWP = new TH1D("hrecotopPt_1m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_1m_MWP->Sumw2();
  hrecotopPt_2m_MWP = new TH1D("hrecotopPt_2m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_2m_MWP->Sumw2();
  hrecotopPt_3m_MWP = new TH1D("hrecotopPt_3m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_3m_MWP->Sumw2();
  hrecotopPt_TWP = new TH1D("hrecotopPt_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_TWP->Sumw2();
  hrecotopPt_0m_TWP = new TH1D("hrecotopPt_0m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_0m_TWP->Sumw2();
  hrecotopPt_1m_TWP = new TH1D("hrecotopPt_1m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_1m_TWP->Sumw2();
  hrecotopPt_2m_TWP = new TH1D("hrecotopPt_2m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_2m_TWP->Sumw2();
  hrecotopPt_3m_TWP = new TH1D("hrecotopPt_3m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_3m_TWP->Sumw2();
  hrecotopPt_XWP = new TH1D("hrecotopPt_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_XWP->Sumw2();
  hrecotopPt_0m_XWP = new TH1D("hrecotopPt_0m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_0m_XWP->Sumw2();
  hrecotopPt_1m_XWP = new TH1D("hrecotopPt_1m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_1m_XWP->Sumw2();
  hrecotopPt_2m_XWP = new TH1D("hrecotopPt_2m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_2m_XWP->Sumw2();
  hrecotopPt_3m_XWP = new TH1D("hrecotopPt_3m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_3m_XWP->Sumw2();
  hrecotopPt = new TH1D("hrecotopPt", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt->Sumw2();
  hrecotopPt_0m = new TH1D("hrecotopPt_0m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_0m->Sumw2();
  hrecotopPt_1m = new TH1D("hrecotopPt_1m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_1m->Sumw2();
  hrecotopPt_2m = new TH1D("hrecotopPt_2m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_2m->Sumw2();
  hrecotopPt_3m = new TH1D("hrecotopPt_3m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_3m->Sumw2();

  hrecotopPt_Mis_LWP = new TH1D("hrecotopPt_Mis_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_LWP->Sumw2();
  hrecotopPt_Mis_0m_LWP = new TH1D("hrecotopPt_Mis_0m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_0m_LWP->Sumw2();
  hrecotopPt_Mis_1m_LWP = new TH1D("hrecotopPt_Mis_1m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_1m_LWP->Sumw2();
  hrecotopPt_Mis_2m_LWP = new TH1D("hrecotopPt_Mis_2m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_2m_LWP->Sumw2();
  hrecotopPt_Mis_3m_LWP = new TH1D("hrecotopPt_Mis_3m_LWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_3m_LWP->Sumw2();
  hrecotopPt_Mis_MWP = new TH1D("hrecotopPt_Mis_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_MWP->Sumw2();
  hrecotopPt_Mis_0m_MWP = new TH1D("hrecotopPt_Mis_0m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_0m_MWP->Sumw2();
  hrecotopPt_Mis_1m_MWP = new TH1D("hrecotopPt_Mis_1m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_1m_MWP->Sumw2();
  hrecotopPt_Mis_2m_MWP = new TH1D("hrecotopPt_Mis_2m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_2m_MWP->Sumw2();
  hrecotopPt_Mis_3m_MWP = new TH1D("hrecotopPt_Mis_3m_MWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_3m_MWP->Sumw2();
  hrecotopPt_Mis_TWP = new TH1D("hrecotopPt_Mis_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_TWP->Sumw2();
  hrecotopPt_Mis_0m_TWP = new TH1D("hrecotopPt_Mis_0m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_0m_TWP->Sumw2();
  hrecotopPt_Mis_1m_TWP = new TH1D("hrecotopPt_Mis_1m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_1m_TWP->Sumw2();
  hrecotopPt_Mis_2m_TWP = new TH1D("hrecotopPt_Mis_2m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_2m_TWP->Sumw2();
  hrecotopPt_Mis_3m_TWP = new TH1D("hrecotopPt_Mis_3m_TWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_3m_TWP->Sumw2();
  hrecotopPt_Mis_XWP = new TH1D("hrecotopPt_Mis_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_XWP->Sumw2();
  hrecotopPt_Mis_0m_XWP = new TH1D("hrecotopPt_Mis_0m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_0m_XWP->Sumw2();
  hrecotopPt_Mis_1m_XWP = new TH1D("hrecotopPt_Mis_1m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_1m_XWP->Sumw2();
  hrecotopPt_Mis_2m_XWP = new TH1D("hrecotopPt_Mis_2m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_2m_XWP->Sumw2();
  hrecotopPt_Mis_3m_XWP = new TH1D("hrecotopPt_Mis_3m_XWP", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_3m_XWP->Sumw2();
  hrecotopPt_Mis = new TH1D("hrecotopPt_Mis", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis->Sumw2();
  hrecotopPt_Mis_0m = new TH1D("hrecotopPt_Mis_0m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_0m->Sumw2();
  hrecotopPt_Mis_1m = new TH1D("hrecotopPt_Mis_1m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_1m->Sumw2();
  hrecotopPt_Mis_2m = new TH1D("hrecotopPt_Mis_2m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_2m->Sumw2();
  hrecotopPt_Mis_3m = new TH1D("hrecotopPt_Mis_3m", "Top p_{T};p_{T}[GeV];Event",20, 0, 1000);
  hrecotopPt_Mis_3m->Sumw2();

  hNjetEff_LWP = new TH1D("hNjetEff_LWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_LWP->Sumw2();
  hNjetEff_0m_LWP = new TH1D("hNjetEff_0m_LWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_0m_LWP->Sumw2();
  hNjetEff_1m_LWP = new TH1D("hNjetEff_1m_LWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_1m_LWP->Sumw2();
  hNjetEff_2m_LWP = new TH1D("hNjetEff_2m_LWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_2m_LWP->Sumw2();
  hNjetEff_3m_LWP = new TH1D("hNjetEff_3m_LWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_3m_LWP->Sumw2();
  hNjetEff_MWP = new TH1D("hNjetEff_MWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_MWP->Sumw2();
  hNjetEff_0m_MWP = new TH1D("hNjetEff_0m_MWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_0m_MWP->Sumw2();
  hNjetEff_1m_MWP = new TH1D("hNjetEff_1m_MWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_1m_MWP->Sumw2();
  hNjetEff_2m_MWP = new TH1D("hNjetEff_2m_MWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_2m_MWP->Sumw2();
  hNjetEff_3m_MWP = new TH1D("hNjetEff_3m_MWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_3m_MWP->Sumw2();
  hNjetEff_TWP = new TH1D("hNjetEff_TWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_TWP->Sumw2();
  hNjetEff_0m_TWP = new TH1D("hNjetEff_0m_TWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_0m_TWP->Sumw2();
  hNjetEff_1m_TWP = new TH1D("hNjetEff_1m_TWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_1m_TWP->Sumw2();
  hNjetEff_2m_TWP = new TH1D("hNjetEff_2m_TWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_2m_TWP->Sumw2();
  hNjetEff_3m_TWP = new TH1D("hNjetEff_3m_TWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_3m_TWP->Sumw2();
  hNjetEff_XWP = new TH1D("hNjetEff_XWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_XWP->Sumw2();
  hNjetEff_0m_XWP = new TH1D("hNjetEff_0m_XWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_0m_XWP->Sumw2();
  hNjetEff_1m_XWP = new TH1D("hNjetEff_1m_XWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_1m_XWP->Sumw2();
  hNjetEff_2m_XWP = new TH1D("hNjetEff_2m_XWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_2m_XWP->Sumw2();
  hNjetEff_3m_XWP = new TH1D("hNjetEff_3m_XWP", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_3m_XWP->Sumw2();
  hNjetEff = new TH1D("hNjetEff", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff->Sumw2();
  hNjetEff_0m = new TH1D("hNjetEff_0m", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_0m->Sumw2();
  hNjetEff_1m = new TH1D("hNjetEff_1m", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_1m->Sumw2();
  hNjetEff_2m = new TH1D("hNjetEff_2m", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_2m->Sumw2();
  hNjetEff_3m = new TH1D("hNjetEff_3m", "N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_3m->Sumw2();
}

bool FillChain(TChain* &chain, const char *subsample, const bool iscondor, const int& startfile, const int& filerun){
  AnaSamples::SampleSet        allSamples = AnaSamples::SampleSet("sampleSets_PostProcessed_2017.cfg", iscondor, AnaSamples::luminosity_muon_2017);
  AnaSamples::SampleCollection allCollections("sampleCollections_2017.cfg", allSamples);
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
