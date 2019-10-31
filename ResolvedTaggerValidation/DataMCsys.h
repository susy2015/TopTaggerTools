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

  TH1D *htopPtEff_den_BTUp;
  TH1D *htopPtEff_den_BTDown;
  TH1D *htopPtEff_den_PUUp;
  TH1D *htopPtEff_den_PUDown;
  TH1D *htopPtEff_den_PDFUp;
  TH1D *htopPtEff_den_PDFDown;
  TH1D *htopPtEff_den_JECUp;
  TH1D *htopPtEff_den_JECDown;
  TH1D *htopPtEff_num_LWP_BTUp;
  TH1D *htopPtEff_num_LWP_BTDown;
  TH1D *htopPtEff_num_LWP_PUUp;
  TH1D *htopPtEff_num_LWP_PUDown;
  TH1D *htopPtEff_num_LWP_PDFUp;
  TH1D *htopPtEff_num_LWP_PDFDown;
  TH1D *htopPtEff_num_LWP_JECUp;
  TH1D *htopPtEff_num_LWP_JECDown;
  TH1D *htopPtEff_num_MWP_BTUp;
  TH1D *htopPtEff_num_MWP_BTDown;
  TH1D *htopPtEff_num_MWP_PUUp;
  TH1D *htopPtEff_num_MWP_PUDown;
  TH1D *htopPtEff_num_MWP_PDFUp;
  TH1D *htopPtEff_num_MWP_PDFDown;
  TH1D *htopPtEff_num_MWP_JECUp;
  TH1D *htopPtEff_num_MWP_JECDown;
  TH1D *htopPtEff_num_TWP_BTUp;
  TH1D *htopPtEff_num_TWP_BTDown;
  TH1D *htopPtEff_num_TWP_PUUp;
  TH1D *htopPtEff_num_TWP_PUDown;
  TH1D *htopPtEff_num_TWP_PDFUp;
  TH1D *htopPtEff_num_TWP_PDFDown;
  TH1D *htopPtEff_num_TWP_JECUp;
  TH1D *htopPtEff_num_TWP_JECDown;
  TH1D *htopPtEff_num_XWP_BTUp;
  TH1D *htopPtEff_num_XWP_BTDown;
  TH1D *htopPtEff_num_XWP_PUUp;
  TH1D *htopPtEff_num_XWP_PUDown;
  TH1D *htopPtEff_num_XWP_PDFUp;
  TH1D *htopPtEff_num_XWP_PDFDown;
  TH1D *htopPtEff_num_XWP_JECUp;
  TH1D *htopPtEff_num_XWP_JECDown;

  TH1D *hNjetEff_den_BTUp;
  TH1D *hNjetEff_den_BTDown;
  TH1D *hNjetEff_den_PUUp;
  TH1D *hNjetEff_den_PUDown;
  TH1D *hNjetEff_den_PDFUp;
  TH1D *hNjetEff_den_PDFDown;
  TH1D *hNjetEff_den_JECUp;
  TH1D *hNjetEff_den_JECDown;
  TH1D *hNjetEff_num_LWP_BTUp;
  TH1D *hNjetEff_num_LWP_BTDown;
  TH1D *hNjetEff_num_LWP_PUUp;
  TH1D *hNjetEff_num_LWP_PUDown;
  TH1D *hNjetEff_num_LWP_PDFUp;
  TH1D *hNjetEff_num_LWP_PDFDown;
  TH1D *hNjetEff_num_LWP_JECUp;
  TH1D *hNjetEff_num_LWP_JECDown;
  TH1D *hNjetEff_num_MWP_BTUp;
  TH1D *hNjetEff_num_MWP_BTDown;
  TH1D *hNjetEff_num_MWP_PUUp;
  TH1D *hNjetEff_num_MWP_PUDown;
  TH1D *hNjetEff_num_MWP_PDFUp;
  TH1D *hNjetEff_num_MWP_PDFDown;
  TH1D *hNjetEff_num_MWP_JECUp;
  TH1D *hNjetEff_num_MWP_JECDown;
  TH1D *hNjetEff_num_TWP_BTUp;
  TH1D *hNjetEff_num_TWP_BTDown;
  TH1D *hNjetEff_num_TWP_PUUp;
  TH1D *hNjetEff_num_TWP_PUDown;
  TH1D *hNjetEff_num_TWP_PDFUp;
  TH1D *hNjetEff_num_TWP_PDFDown;
  TH1D *hNjetEff_num_TWP_JECUp;
  TH1D *hNjetEff_num_TWP_JECDown;
  TH1D *hNjetEff_num_XWP_BTUp;
  TH1D *hNjetEff_num_XWP_BTDown;
  TH1D *hNjetEff_num_XWP_PUUp;
  TH1D *hNjetEff_num_XWP_PUDown;
  TH1D *hNjetEff_num_XWP_PDFUp;
  TH1D *hNjetEff_num_XWP_PDFDown;
  TH1D *hNjetEff_num_XWP_JECUp;
  TH1D *hNjetEff_num_XWP_JECDown;

  TH1D *htopPtMis_den_BTUp;
  TH1D *htopPtMis_den_BTDown;
  TH1D *htopPtMis_den_PUUp;
  TH1D *htopPtMis_den_PUDown;
  TH1D *htopPtMis_den_PDFUp;
  TH1D *htopPtMis_den_PDFDown;
  TH1D *htopPtMis_den_JECUp;
  TH1D *htopPtMis_den_JECDown;
  TH1D *htopPtMis_num_LWP_BTUp;
  TH1D *htopPtMis_num_LWP_BTDown;
  TH1D *htopPtMis_num_LWP_PUUp;
  TH1D *htopPtMis_num_LWP_PUDown;
  TH1D *htopPtMis_num_LWP_PDFUp;
  TH1D *htopPtMis_num_LWP_PDFDown;
  TH1D *htopPtMis_num_LWP_JECUp;
  TH1D *htopPtMis_num_LWP_JECDown;
  TH1D *htopPtMis_num_MWP_BTUp;
  TH1D *htopPtMis_num_MWP_BTDown;
  TH1D *htopPtMis_num_MWP_PUUp;
  TH1D *htopPtMis_num_MWP_PUDown;
  TH1D *htopPtMis_num_MWP_PDFUp;
  TH1D *htopPtMis_num_MWP_PDFDown;
  TH1D *htopPtMis_num_MWP_JECUp;
  TH1D *htopPtMis_num_MWP_JECDown;
  TH1D *htopPtMis_num_TWP_BTUp;
  TH1D *htopPtMis_num_TWP_BTDown;
  TH1D *htopPtMis_num_TWP_PUUp;
  TH1D *htopPtMis_num_TWP_PUDown;
  TH1D *htopPtMis_num_TWP_PDFUp;
  TH1D *htopPtMis_num_TWP_PDFDown;
  TH1D *htopPtMis_num_TWP_JECUp;
  TH1D *htopPtMis_num_TWP_JECDown;
  TH1D *htopPtMis_num_XWP_BTUp;
  TH1D *htopPtMis_num_XWP_BTDown;
  TH1D *htopPtMis_num_XWP_PUUp;
  TH1D *htopPtMis_num_XWP_PUDown;
  TH1D *htopPtMis_num_XWP_PDFUp;
  TH1D *htopPtMis_num_XWP_PDFDown;
  TH1D *htopPtMis_num_XWP_JECUp;
  TH1D *htopPtMis_num_XWP_JECDown;

  TH1D *hNjetMis_den_BTUp;
  TH1D *hNjetMis_den_BTDown;
  TH1D *hNjetMis_den_PUUp;
  TH1D *hNjetMis_den_PUDown;
  TH1D *hNjetMis_den_PDFUp;
  TH1D *hNjetMis_den_PDFDown;
  TH1D *hNjetMis_den_JECUp;
  TH1D *hNjetMis_den_JECDown;
  TH1D *hNjetMis_num_LWP_BTUp;
  TH1D *hNjetMis_num_LWP_BTDown;
  TH1D *hNjetMis_num_LWP_PUUp;
  TH1D *hNjetMis_num_LWP_PUDown;
  TH1D *hNjetMis_num_LWP_PDFUp;
  TH1D *hNjetMis_num_LWP_PDFDown;
  TH1D *hNjetMis_num_LWP_JECUp;
  TH1D *hNjetMis_num_LWP_JECDown;
  TH1D *hNjetMis_num_MWP_BTUp;
  TH1D *hNjetMis_num_MWP_BTDown;
  TH1D *hNjetMis_num_MWP_PUUp;
  TH1D *hNjetMis_num_MWP_PUDown;
  TH1D *hNjetMis_num_MWP_PDFUp;
  TH1D *hNjetMis_num_MWP_PDFDown;
  TH1D *hNjetMis_num_MWP_JECUp;
  TH1D *hNjetMis_num_MWP_JECDown;
  TH1D *hNjetMis_num_TWP_BTUp;
  TH1D *hNjetMis_num_TWP_BTDown;
  TH1D *hNjetMis_num_TWP_PUUp;
  TH1D *hNjetMis_num_TWP_PUDown;
  TH1D *hNjetMis_num_TWP_PDFUp;
  TH1D *hNjetMis_num_TWP_PDFDown;
  TH1D *hNjetMis_num_TWP_JECUp;
  TH1D *hNjetMis_num_TWP_JECDown;
  TH1D *hNjetMis_num_XWP_BTUp;
  TH1D *hNjetMis_num_XWP_BTDown;
  TH1D *hNjetMis_num_XWP_PUUp;
  TH1D *hNjetMis_num_XWP_PUDown;
  TH1D *hNjetMis_num_XWP_PDFUp;
  TH1D *hNjetMis_num_XWP_PDFDown;
  TH1D *hNjetMis_num_XWP_JECUp;
  TH1D *hNjetMis_num_XWP_JECDown;
  
};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_DataMCsys"+index+".root";
  oFile = new TFile(filename, "recreate");
 
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

  htopPtEff_den_BTUp = new TH1D("htopPtEff_den_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_BTUp->Sumw2();
  htopPtEff_den_BTDown = new TH1D("htopPtEff_den_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_BTDown->Sumw2();
  htopPtEff_den_PUUp = new TH1D("htopPtEff_den_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_PUUp->Sumw2();
  htopPtEff_den_PUDown = new TH1D("htopPtEff_den_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_PUDown->Sumw2();
  htopPtEff_den_PDFUp = new TH1D("htopPtEff_den_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_PDFUp->Sumw2();
  htopPtEff_den_PDFDown = new TH1D("htopPtEff_den_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_PDFDown->Sumw2();
  htopPtEff_den_JECUp = new TH1D("htopPtEff_den_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_JECUp->Sumw2();
  htopPtEff_den_JECDown = new TH1D("htopPtEff_den_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_den_JECDown->Sumw2();
  htopPtEff_num_LWP_BTUp = new TH1D("htopPtEff_num_LWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_BTUp->Sumw2();
  htopPtEff_num_LWP_BTDown = new TH1D("htopPtEff_num_LWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_BTDown->Sumw2();
  htopPtEff_num_LWP_PUUp = new TH1D("htopPtEff_num_LWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_PUUp->Sumw2();
  htopPtEff_num_LWP_PUDown = new TH1D("htopPtEff_num_LWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_PUDown->Sumw2();
  htopPtEff_num_LWP_PDFUp = new TH1D("htopPtEff_num_LWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_PDFUp->Sumw2();
  htopPtEff_num_LWP_PDFDown = new TH1D("htopPtEff_num_LWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_PDFDown->Sumw2();
  htopPtEff_num_LWP_JECUp = new TH1D("htopPtEff_num_LWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_JECUp->Sumw2();
  htopPtEff_num_LWP_JECDown = new TH1D("htopPtEff_num_LWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_LWP_JECDown->Sumw2();
  htopPtEff_num_MWP_BTUp = new TH1D("htopPtEff_num_MWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_BTUp->Sumw2();
  htopPtEff_num_MWP_BTDown = new TH1D("htopPtEff_num_MWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_BTDown->Sumw2();
  htopPtEff_num_MWP_PUUp = new TH1D("htopPtEff_num_MWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_PUUp->Sumw2();
  htopPtEff_num_MWP_PUDown = new TH1D("htopPtEff_num_MWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_PUDown->Sumw2();
  htopPtEff_num_MWP_PDFUp = new TH1D("htopPtEff_num_MWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_PDFUp->Sumw2();
  htopPtEff_num_MWP_PDFDown = new TH1D("htopPtEff_num_MWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_PDFDown->Sumw2();
  htopPtEff_num_MWP_JECUp = new TH1D("htopPtEff_num_MWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_JECUp->Sumw2();
  htopPtEff_num_MWP_JECDown = new TH1D("htopPtEff_num_MWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_MWP_JECDown->Sumw2();
  htopPtEff_num_TWP_BTUp = new TH1D("htopPtEff_num_TWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_BTUp->Sumw2();
  htopPtEff_num_TWP_BTDown = new TH1D("htopPtEff_num_TWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_BTDown->Sumw2();
  htopPtEff_num_TWP_PUUp = new TH1D("htopPtEff_num_TWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_PUUp->Sumw2();
  htopPtEff_num_TWP_PUDown = new TH1D("htopPtEff_num_TWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_PUDown->Sumw2();
  htopPtEff_num_TWP_PDFUp = new TH1D("htopPtEff_num_TWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_PDFUp->Sumw2();
  htopPtEff_num_TWP_PDFDown = new TH1D("htopPtEff_num_TWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_PDFDown->Sumw2();
  htopPtEff_num_TWP_JECUp = new TH1D("htopPtEff_num_TWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_JECUp->Sumw2();
  htopPtEff_num_TWP_JECDown = new TH1D("htopPtEff_num_TWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_TWP_JECDown->Sumw2();
  htopPtEff_num_XWP_BTUp = new TH1D("htopPtEff_num_XWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_BTUp->Sumw2();
  htopPtEff_num_XWP_BTDown = new TH1D("htopPtEff_num_XWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_BTDown->Sumw2();
  htopPtEff_num_XWP_PUUp = new TH1D("htopPtEff_num_XWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_PUUp->Sumw2();
  htopPtEff_num_XWP_PUDown = new TH1D("htopPtEff_num_XWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_PUDown->Sumw2();
  htopPtEff_num_XWP_PDFUp = new TH1D("htopPtEff_num_XWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_PDFUp->Sumw2();
  htopPtEff_num_XWP_PDFDown = new TH1D("htopPtEff_num_XWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_PDFDown->Sumw2();
  htopPtEff_num_XWP_JECUp = new TH1D("htopPtEff_num_XWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_JECUp->Sumw2();
  htopPtEff_num_XWP_JECDown = new TH1D("htopPtEff_num_XWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtEff_num_XWP_JECDown->Sumw2();

  hNjetEff_den_BTUp = new TH1D("hNjetEff_den_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_BTUp->Sumw2();
  hNjetEff_den_BTDown = new TH1D("hNjetEff_den_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_BTDown->Sumw2();
  hNjetEff_den_PUUp = new TH1D("hNjetEff_den_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_PUUp->Sumw2();
  hNjetEff_den_PUDown = new TH1D("hNjetEff_den_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_PUDown->Sumw2();
  hNjetEff_den_PDFUp = new TH1D("hNjetEff_den_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_PDFUp->Sumw2();
  hNjetEff_den_PDFDown = new TH1D("hNjetEff_den_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_PDFDown->Sumw2();
  hNjetEff_den_JECUp = new TH1D("hNjetEff_den_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_JECUp->Sumw2();
  hNjetEff_den_JECDown = new TH1D("hNjetEff_den_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_den_JECDown->Sumw2();
  hNjetEff_num_LWP_BTUp = new TH1D("hNjetEff_num_LWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_BTUp->Sumw2();
  hNjetEff_num_LWP_BTDown = new TH1D("hNjetEff_num_LWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_BTDown->Sumw2();
  hNjetEff_num_LWP_PUUp = new TH1D("hNjetEff_num_LWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_PUUp->Sumw2();
  hNjetEff_num_LWP_PUDown = new TH1D("hNjetEff_num_LWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_PUDown->Sumw2();
  hNjetEff_num_LWP_PDFUp = new TH1D("hNjetEff_num_LWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_PDFUp->Sumw2();
  hNjetEff_num_LWP_PDFDown = new TH1D("hNjetEff_num_LWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_PDFDown->Sumw2();
  hNjetEff_num_LWP_JECUp = new TH1D("hNjetEff_num_LWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_JECUp->Sumw2();
  hNjetEff_num_LWP_JECDown = new TH1D("hNjetEff_num_LWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_LWP_JECDown->Sumw2();
  hNjetEff_num_MWP_BTUp = new TH1D("hNjetEff_num_MWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_BTUp->Sumw2();
  hNjetEff_num_MWP_BTDown = new TH1D("hNjetEff_num_MWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_BTDown->Sumw2();
  hNjetEff_num_MWP_PUUp = new TH1D("hNjetEff_num_MWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_PUUp->Sumw2();
  hNjetEff_num_MWP_PUDown = new TH1D("hNjetEff_num_MWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_PUDown->Sumw2();
  hNjetEff_num_MWP_PDFUp = new TH1D("hNjetEff_num_MWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_PDFUp->Sumw2();
  hNjetEff_num_MWP_PDFDown = new TH1D("hNjetEff_num_MWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_PDFDown->Sumw2();
  hNjetEff_num_MWP_JECUp = new TH1D("hNjetEff_num_MWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_JECUp->Sumw2();
  hNjetEff_num_MWP_JECDown = new TH1D("hNjetEff_num_MWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_MWP_JECDown->Sumw2();
  hNjetEff_num_TWP_BTUp = new TH1D("hNjetEff_num_TWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_BTUp->Sumw2();
  hNjetEff_num_TWP_BTDown = new TH1D("hNjetEff_num_TWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_BTDown->Sumw2();
  hNjetEff_num_TWP_PUUp = new TH1D("hNjetEff_num_TWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_PUUp->Sumw2();
  hNjetEff_num_TWP_PUDown = new TH1D("hNjetEff_num_TWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_PUDown->Sumw2();
  hNjetEff_num_TWP_PDFUp = new TH1D("hNjetEff_num_TWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_PDFUp->Sumw2();
  hNjetEff_num_TWP_PDFDown = new TH1D("hNjetEff_num_TWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_PDFDown->Sumw2();
  hNjetEff_num_TWP_JECUp = new TH1D("hNjetEff_num_TWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_JECUp->Sumw2();
  hNjetEff_num_TWP_JECDown = new TH1D("hNjetEff_num_TWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_TWP_JECDown->Sumw2();
  hNjetEff_num_XWP_BTUp = new TH1D("hNjetEff_num_XWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_BTUp->Sumw2();
  hNjetEff_num_XWP_BTDown = new TH1D("hNjetEff_num_XWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_BTDown->Sumw2();
  hNjetEff_num_XWP_PUUp = new TH1D("hNjetEff_num_XWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_PUUp->Sumw2();
  hNjetEff_num_XWP_PUDown = new TH1D("hNjetEff_num_XWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_PUDown->Sumw2();
  hNjetEff_num_XWP_PDFUp = new TH1D("hNjetEff_num_XWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_PDFUp->Sumw2();
  hNjetEff_num_XWP_PDFDown = new TH1D("hNjetEff_num_XWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_PDFDown->Sumw2();
  hNjetEff_num_XWP_JECUp = new TH1D("hNjetEff_num_XWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_JECUp->Sumw2();
  hNjetEff_num_XWP_JECDown = new TH1D("hNjetEff_num_XWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetEff_num_XWP_JECDown->Sumw2();

  htopPtMis_den_BTUp = new TH1D("htopPtMis_den_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_BTUp->Sumw2();
  htopPtMis_den_BTDown = new TH1D("htopPtMis_den_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_BTDown->Sumw2();
  htopPtMis_den_PUUp = new TH1D("htopPtMis_den_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_PUUp->Sumw2();
  htopPtMis_den_PUDown = new TH1D("htopPtMis_den_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_PUDown->Sumw2();
  htopPtMis_den_PDFUp = new TH1D("htopPtMis_den_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_PDFUp->Sumw2();
  htopPtMis_den_PDFDown = new TH1D("htopPtMis_den_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_PDFDown->Sumw2();
  htopPtMis_den_JECUp = new TH1D("htopPtMis_den_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_JECUp->Sumw2();
  htopPtMis_den_JECDown = new TH1D("htopPtMis_den_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_den_JECDown->Sumw2();
  htopPtMis_num_LWP_BTUp = new TH1D("htopPtMis_num_LWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_BTUp->Sumw2();
  htopPtMis_num_LWP_BTDown = new TH1D("htopPtMis_num_LWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_BTDown->Sumw2();
  htopPtMis_num_LWP_PUUp = new TH1D("htopPtMis_num_LWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_PUUp->Sumw2();
  htopPtMis_num_LWP_PUDown = new TH1D("htopPtMis_num_LWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_PUDown->Sumw2();
  htopPtMis_num_LWP_PDFUp = new TH1D("htopPtMis_num_LWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_PDFUp->Sumw2();
  htopPtMis_num_LWP_PDFDown = new TH1D("htopPtMis_num_LWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_PDFDown->Sumw2();
  htopPtMis_num_LWP_JECUp = new TH1D("htopPtMis_num_LWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_JECUp->Sumw2();
  htopPtMis_num_LWP_JECDown = new TH1D("htopPtMis_num_LWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_LWP_JECDown->Sumw2();
  htopPtMis_num_MWP_BTUp = new TH1D("htopPtMis_num_MWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_BTUp->Sumw2();
  htopPtMis_num_MWP_BTDown = new TH1D("htopPtMis_num_MWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_BTDown->Sumw2();
  htopPtMis_num_MWP_PUUp = new TH1D("htopPtMis_num_MWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_PUUp->Sumw2();
  htopPtMis_num_MWP_PUDown = new TH1D("htopPtMis_num_MWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_PUDown->Sumw2();
  htopPtMis_num_MWP_PDFUp = new TH1D("htopPtMis_num_MWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_PDFUp->Sumw2();
  htopPtMis_num_MWP_PDFDown = new TH1D("htopPtMis_num_MWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_PDFDown->Sumw2();
  htopPtMis_num_MWP_JECUp = new TH1D("htopPtMis_num_MWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_JECUp->Sumw2();
  htopPtMis_num_MWP_JECDown = new TH1D("htopPtMis_num_MWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_MWP_JECDown->Sumw2();
  htopPtMis_num_TWP_BTUp = new TH1D("htopPtMis_num_TWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_BTUp->Sumw2();
  htopPtMis_num_TWP_BTDown = new TH1D("htopPtMis_num_TWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_BTDown->Sumw2();
  htopPtMis_num_TWP_PUUp = new TH1D("htopPtMis_num_TWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_PUUp->Sumw2();
  htopPtMis_num_TWP_PUDown = new TH1D("htopPtMis_num_TWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_PUDown->Sumw2();
  htopPtMis_num_TWP_PDFUp = new TH1D("htopPtMis_num_TWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_PDFUp->Sumw2();
  htopPtMis_num_TWP_PDFDown = new TH1D("htopPtMis_num_TWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_PDFDown->Sumw2();
  htopPtMis_num_TWP_JECUp = new TH1D("htopPtMis_num_TWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_JECUp->Sumw2();
  htopPtMis_num_TWP_JECDown = new TH1D("htopPtMis_num_TWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_TWP_JECDown->Sumw2();
  htopPtMis_num_XWP_BTUp = new TH1D("htopPtMis_num_XWP_BTUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_BTUp->Sumw2();
  htopPtMis_num_XWP_BTDown = new TH1D("htopPtMis_num_XWP_BTDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_BTDown->Sumw2();
  htopPtMis_num_XWP_PUUp = new TH1D("htopPtMis_num_XWP_PUUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_PUUp->Sumw2();
  htopPtMis_num_XWP_PUDown = new TH1D("htopPtMis_num_XWP_PUDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_PUDown->Sumw2();
  htopPtMis_num_XWP_PDFUp = new TH1D("htopPtMis_num_XWP_PDFUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_PDFUp->Sumw2();
  htopPtMis_num_XWP_PDFDown = new TH1D("htopPtMis_num_XWP_PDFDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_PDFDown->Sumw2();
  htopPtMis_num_XWP_JECUp = new TH1D("htopPtMis_num_XWP_JECUp", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_JECUp->Sumw2();
  htopPtMis_num_XWP_JECDown = new TH1D("htopPtMis_num_XWP_JECDown", "Top (had) p_{T};p_{T}[GeV];Event",20, 0, 1000);
  htopPtMis_num_XWP_JECDown->Sumw2();

  hNjetMis_den_BTUp = new TH1D("hNjetMis_den_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_BTUp->Sumw2();
  hNjetMis_den_BTDown = new TH1D("hNjetMis_den_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_BTDown->Sumw2();
  hNjetMis_den_PUUp = new TH1D("hNjetMis_den_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_PUUp->Sumw2();
  hNjetMis_den_PUDown = new TH1D("hNjetMis_den_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_PUDown->Sumw2();
  hNjetMis_den_PDFUp = new TH1D("hNjetMis_den_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_PDFUp->Sumw2();
  hNjetMis_den_PDFDown = new TH1D("hNjetMis_den_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_PDFDown->Sumw2();  
  hNjetMis_den_JECUp = new TH1D("hNjetMis_den_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_JECUp->Sumw2();
  hNjetMis_den_JECDown = new TH1D("hNjetMis_den_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_den_JECDown->Sumw2();
  
  
  hNjetMis_num_LWP_BTUp = new TH1D("hNjetMis_num_LWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_BTUp->Sumw2();
  hNjetMis_num_LWP_BTDown = new TH1D("hNjetMis_num_LWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_BTDown->Sumw2();
  hNjetMis_num_LWP_PUUp = new TH1D("hNjetMis_num_LWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_PUUp->Sumw2();
  hNjetMis_num_LWP_PUDown = new TH1D("hNjetMis_num_LWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_PUDown->Sumw2();
  hNjetMis_num_LWP_PDFUp = new TH1D("hNjetMis_num_LWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_PDFUp->Sumw2();
  hNjetMis_num_LWP_PDFDown = new TH1D("hNjetMis_num_LWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_PDFDown->Sumw2();
  hNjetMis_num_LWP_JECUp = new TH1D("hNjetMis_num_LWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_JECUp->Sumw2();
  hNjetMis_num_LWP_JECDown = new TH1D("hNjetMis_num_LWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_LWP_JECDown->Sumw2();
  hNjetMis_num_MWP_BTUp = new TH1D("hNjetMis_num_MWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_BTUp->Sumw2();
  hNjetMis_num_MWP_BTDown = new TH1D("hNjetMis_num_MWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_BTDown->Sumw2();
  hNjetMis_num_MWP_PUUp = new TH1D("hNjetMis_num_MWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_PUUp->Sumw2();
  hNjetMis_num_MWP_PUDown = new TH1D("hNjetMis_num_MWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_PUDown->Sumw2();
  hNjetMis_num_MWP_PDFUp = new TH1D("hNjetMis_num_MWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_PDFUp->Sumw2();
  hNjetMis_num_MWP_PDFDown = new TH1D("hNjetMis_num_MWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_PDFDown->Sumw2();
  hNjetMis_num_MWP_JECUp = new TH1D("hNjetMis_num_MWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_JECUp->Sumw2();
  hNjetMis_num_MWP_JECDown = new TH1D("hNjetMis_num_MWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_MWP_JECDown->Sumw2();
  hNjetMis_num_TWP_BTUp = new TH1D("hNjetMis_num_TWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_BTUp->Sumw2();
  hNjetMis_num_TWP_BTDown = new TH1D("hNjetMis_num_TWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_BTDown->Sumw2();
  hNjetMis_num_TWP_PUUp = new TH1D("hNjetMis_num_TWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_PUUp->Sumw2();
  hNjetMis_num_TWP_PUDown = new TH1D("hNjetMis_num_TWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_PUDown->Sumw2();
  hNjetMis_num_TWP_PDFUp = new TH1D("hNjetMis_num_TWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_PDFUp->Sumw2();
  hNjetMis_num_TWP_PDFDown = new TH1D("hNjetMis_num_TWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_PDFDown->Sumw2();
  hNjetMis_num_TWP_JECUp = new TH1D("hNjetMis_num_TWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_JECUp->Sumw2();
  hNjetMis_num_TWP_JECDown = new TH1D("hNjetMis_num_TWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_TWP_JECDown->Sumw2();
  hNjetMis_num_XWP_BTUp = new TH1D("hNjetMis_num_XWP_BTUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_BTUp->Sumw2();
  hNjetMis_num_XWP_BTDown = new TH1D("hNjetMis_num_XWP_BTDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_BTDown->Sumw2();
  hNjetMis_num_XWP_PUUp = new TH1D("hNjetMis_num_XWP_PUUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_PUUp->Sumw2();
  hNjetMis_num_XWP_PUDown = new TH1D("hNjetMis_num_XWP_PUDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_PUDown->Sumw2();
  hNjetMis_num_XWP_PDFUp = new TH1D("hNjetMis_num_XWP_PDFUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_PDFUp->Sumw2();
  hNjetMis_num_XWP_PDFDown = new TH1D("hNjetMis_num_XWP_PDFDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_PDFDown->Sumw2();
  hNjetMis_num_XWP_JECUp = new TH1D("hNjetMis_num_XWP_JECUp","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_JECUp->Sumw2();
  hNjetMis_num_XWP_JECDown = new TH1D("hNjetMis_num_XWP_JECDown","N_{jets};N_{jet};Event",15,0,15);
  hNjetMis_num_XWP_JECDown->Sumw2();

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
