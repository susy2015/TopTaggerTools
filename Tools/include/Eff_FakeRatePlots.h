#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include <TLegend.h>
#include <TCanvas.h>
#include "TF1.h"

#include <iostream>
#include <string>

#include <cstdlib>
#include <cstdio>
#include <cstring>

class Eff_FakeRatePlots
{

private:  
  
public:
    std::vector<TH1*> histos_;
    std::vector<std::string> histoName_;
    std::vector<std::string> outFile_;

    TH1F* defineTH1F(const char *file="", std::string branch="")
    {
        //TFile *f = new TFile(file);
        TFile *f = TFile::Open(file);
        TH1F *h = (TH1F*)f->Get(branch.c_str());
        f->Close();
        return h;
    }

    double errorRatio(int j, TH1* hnum, TH1* hdom)
    {
        double a     = hnum->GetBinContent(j);
        double b     = hdom->GetBinContent(j);
        double da    = hnum->GetBinError(j);
        double db    = hdom->GetBinError(j);    
        double error = (a/b)*sqrt( pow(da/a,2) + pow(db/b,2) );     
        if(isinf(error)==1 || isnan(error)==1) error = 0;    
        return error;
    }
  
    void ratioTH1F(std::string name, std::string fakerateVar, std::string Var, const char* file, std::string fakerate, std::string type, std::string xname, std::string yname, int join )
    {
        TH1F* fakeRateVar = defineTH1F(file,fakerateVar);
        TH1F* var         = defineTH1F(file,Var);
        fakeRateVar->Rebin(join);
        var->Rebin(join);

        TH1F* fakeRate    = (TH1F*)fakeRateVar->Clone(fakerate.c_str());
        for(int i=0; i < fakeRate->GetSize(); i++)
        {
            double error = errorRatio(i,fakeRate,var);
            double ratio = fakeRate->GetBinContent(i)/var->GetBinContent(i);      
            if(isinf(ratio)==1 || isnan(ratio)==1) ratio = 0;
            fakeRate->SetBinContent(i,ratio);
            fakeRate->SetBinError(i,error);
        }
        fakeRate->SetTitle((type + "/" + fakerate).c_str());
        fakeRate->GetXaxis()->SetTitle(xname.c_str());
        fakeRate->GetYaxis()->SetTitle(yname.c_str());
        histos_.emplace_back(fakeRate);
        histoName_.emplace_back(type + "/" + fakerate);
        outFile_.emplace_back(type + "/" + name);
    }

    void makeTH1F(std::string name, const char* rootFile, std::string type)
    {
        ratioTH1F(name,  name + "fakerateMET"     ,name + "MET"               ,rootFile ,name + "fakerateMET"              ,type ,"MET (GeV)"     ,"Fakerate"  , 5);
        ratioTH1F(name,  name + "fakerateNj"      ,name + "nJets"             ,rootFile ,name + "fakerateNjets"            ,type ,"N_{Jets}"      ,"Fakerate"  , 1);
        ratioTH1F(name,  name + "fakerateNb"      ,name + "nBJets"            ,rootFile ,name + "fakerateNbjets"           ,type ,"N_{BJets}"     ,"Fakerate"  , 1);			     
        ratioTH1F(name,  name + "fakerateMET2"    ,name + "MET"               ,rootFile ,name + "fakerateMET2"             ,type ,"MET (GeV)"     ,"Fakerate"  , 5);
        ratioTH1F(name,  name + "fakerateNj2"     ,name + "nJets"             ,rootFile ,name + "fakerateNjets2"           ,type ,"N_{Jets}"      ,"Fakerate"  , 1);
        ratioTH1F(name,  name + "fakerateNb2"     ,name + "nBJets"            ,rootFile ,name + "fakerateNbjets2"          ,type ,"N_{BJets}"     ,"Fakerate"  , 1);    			     
        ratioTH1F(name,  name + "randomTopPt"     ,name + "randomTopCandPt"   ,rootFile ,name + "fakerateRandomTopPt"      ,type ,"P_{T} (GeV)"   ,"Fakerate"  , 5);
        ratioTH1F(name,  name + "randomTopMass"   ,name + "randomTopCandMass" ,rootFile ,name + "fakerateRandomTopMass"    ,type ,"Mass (GeV)"    ,"Fakerate"  , 1);
        ratioTH1F(name,  name + "randomTopEta"    ,name + "randomTopCandEta"  ,rootFile ,name + "fakerateRandomTopEta"     ,type ,"#eta"          ,"Fakerate"  , 2);    
        ratioTH1F(name,  name + "genTopMatchPt"   ,name + "genTopPt"          ,rootFile ,name + "efficiencyGenTopMatchPt"  ,type ,"P_{T} (GeV)"   ,"Efficiency", 5);
        ratioTH1F(name,  name + "genTopMatchMass" ,name + "genTopMass"        ,rootFile ,name + "efficiencyGenTopMatchMass",type ,"Mass (GeV)"    ,"Efficiency", 1);	
        ratioTH1F(name,  name + "genTopMatchEta"  ,name + "genTopEta"         ,rootFile ,name + "efficiencyGenTopMatchEta" ,type ,"#eta"          ,"Efficiency", 2);
    }
    
    Eff_FakeRatePlots(){}
    ~Eff_FakeRatePlots(){}

};
