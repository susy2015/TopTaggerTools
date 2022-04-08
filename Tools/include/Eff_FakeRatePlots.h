#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include <TLegend.h>
#include <TCanvas.h>
#include "TLatex.h"
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
        TFile *f = TFile::Open(file);
        if(f)
        {
            f->cd();
            TH1F *h = (TH1F*)f->Get(branch.c_str());
            f->Close();
            return h;
        }
        else 
        {
            std::cout<<"File is nullptr"<<std::endl;
            return nullptr;
        }
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
        //fakeRate->SetTitle((type + "/" + fakerate).c_str());
        //fakeRate->SetTitle( "Resolved Top Tagger" );
        fakeRate->GetXaxis()->SetTitle(xname.c_str());
        fakeRate->GetYaxis()->SetTitle(yname.c_str());
        fakeRate->GetYaxis()->SetRangeUser(0.0, 1.2); 
        fakeRate->SetStats(0);
        histos_.emplace_back(fakeRate);
        histoName_.emplace_back(type + "/" + fakerate);
        outFile_.emplace_back(type + "/" + name + "/"); 
    }

    void makeTH1F(std::string name, const char* rootFile, std::string type, std::string wp, std::string dataset)
    {
        if (dataset == "TT")
        {
            ratioTH1F(name, name + "/genTopMatchPtWP_" + wp + name, name + "/genTopPt_" + name, rootFile, name + "/efficiencyGenTopMatchPt", type ,"p_{T} (GeV)", "Efficiency", 5);
        } else
        {
            ratioTH1F(name, name + "/fakerateNj_"      + wp + name, name + "/nJets_"    + name, rootFile, name + "/fakerateNjets"          , type ,"N_{jets}"   , "Fakerate"  , 1);
        } 
    }
 
    Eff_FakeRatePlots(){}
    ~Eff_FakeRatePlots(){}
};
