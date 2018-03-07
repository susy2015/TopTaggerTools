#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include <TLegend.h>
#include <TCanvas.h>

#include <iostream>
#include <string>

class FakeRatePlots
{

private:  
  
public:
    std::vector<TH1*> histos_;
    std::vector<std::string> histoName_;
    
    TH1F* defineTH1F(const char *file="", std::string branch="")
    {
        TFile *f = new TFile(file);
        TH1F *h = (TH1F*)f->Get(branch.c_str());
        return h;
    }

    double errorRatio(int j, TH1F* hnum, TH1F* hdom)
    {
        double a     = hnum->GetBinContent(j);
        double b     = hdom->GetBinContent(j);
        double da    = hnum->GetBinError(j);
        double db    = hdom->GetBinError(j);    
        double error = (a/b)*sqrt( pow(da/a,2) + pow(db/b,2) );     
        if(isinf(error)==1 || isnan(error)==1) error = 0;    
        return error;
    }
  
    void ratioTH1F(std::string fakerateVar, std::string Var, const char* file, std::string fakerate)
    {
        TH1F* fakeRateVar = defineTH1F(file,fakerateVar);
        TH1F* var         = defineTH1F(file,Var);
        TH1F* fakeRate    = (TH1F*)fakeRateVar->Clone(fakerate.c_str());
        for(int i=0; i < fakeRate->GetSize(); i++)
        {
            double error = errorRatio(i,fakeRate,var);
            double ratio = fakeRate->GetBinContent(i)/var->GetBinContent(i);      
            if(isinf(ratio)==1 || isnan(ratio)==1) ratio = 0;
            fakeRate->SetBinContent(i,ratio);
            fakeRate->SetBinError(i,error);
        }
        fakeRate->SetTitle(fakerate.c_str());
        histos_.emplace_back(fakeRate);
        histoName_.emplace_back(fakerate);
    }

    void makeTH1F(std::string name, const char* rootFile)
    {
        ratioTH1F( name + "fakerateMET"     ,name + "MET"              ,rootFile, name + "fakerateMET"              );
        ratioTH1F( name + "fakerateNj"      ,name + "nJets"            ,rootFile, name + "fakerateNjets"            );
        ratioTH1F( name + "fakerateNb"      ,name + "nBJets"           ,rootFile, name + "fakerateNbjets"           );			     
        ratioTH1F( name + "fakerateMET2"    ,name + "MET"              ,rootFile, name + "fakerateMET2"             );
        ratioTH1F( name + "fakerateNj2"     ,name + "nJets"            ,rootFile, name + "fakerateNjets2"           );
        ratioTH1F( name + "fakerateNb2"     ,name + "nBJets"           ,rootFile, name + "fakerateNbjets2"          );    			     
        ratioTH1F( name + "randomTopPt"     ,name + "randomTopCandPt"  ,rootFile, name + "fakerateRandomTopPt"      );
        ratioTH1F( name + "randomTopMass"   ,name + "randomTopCandMass",rootFile, name + "fakerateRandomTopMass"    );
        ratioTH1F( name + "randomTopEta"    ,name + "randomTopCandEta" ,rootFile, name + "fakerateRandomTopEta"     );    
        ratioTH1F( name + "genTopMatchPt"   ,name + "genTopPt"         ,rootFile, name + "efficiencyGenTopMatchPt"  );
        ratioTH1F( name + "genTopMatchMass" ,name + "genTopMass"       ,rootFile, name + "efficiencyGenTopMatchMass");	
        ratioTH1F( name + "genTopMatchEta"  ,name + "genTopEta"        ,rootFile, name + "efficiencyGenTopMatchEta" );
    }
    
    FakeRatePlots(){}
    ~FakeRatePlots(){}
  
};


void makePlots(std::string name, TH1* simpleHist, TH1* mediumHist)
{
    TCanvas *c = new TCanvas(name.c_str(),name.c_str(),1000,800);  
    TLegend *l = new TLegend(0.68,0.8,0.99,0.9);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.16);
    
    //mediumHist->GetXaxis()->SetTitle("time (ns)");
    //mediumHist->GetYaxis()->SetTitle("A.U.");
    //mediumHist->SetTitle("Pulse Shape Comparison");
    //mediumHist->SetName("Pulse Shape  Comparison");
    mediumHist->SetTitleSize(0.002);
    mediumHist->SetTitleSize(0.05,"X");
    mediumHist->SetTitleSize(0.04,"Y");
    mediumHist->SetTitleOffset(1.2,"X");
    mediumHist->SetTitleOffset(1.8,"Y");
    mediumHist->SetLabelSize(0.05,"X");
    mediumHist->SetLabelSize(0.05,"Y");
    mediumHist->SetStats(false);    
    mediumHist->SetLineColor(kRed);
    mediumHist->Draw("hist E");

    simpleHist->SetLineColor(kBlack);
    simpleHist->Draw("hist E same");

    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(0.03);    
    l->AddEntry(mediumHist, "Medium Top Tagger", "l");
    l->AddEntry(simpleHist, "Simple   Top Tagger", "l");
    l->Draw();
    
    c->SaveAs(("plots/" + name + ".png").c_str());        
}

int main()
{    
    ///////////////////////////////////
    //Simple Top Tagger
    ///////////////////////////////////  
    const char* rootFileSimple = "joesGroup/TT_TTbarSingleLep_simpleTopTagger.root";
    std::string filenameSimple = "efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root";
    
    FakeRatePlots fakeratePlotsSimple;    
    fakeratePlotsSimple.makeTH1F("Lep0_", rootFileSimple);
    fakeratePlotsSimple.makeTH1F("Lep1_", rootFileSimple);

    TFile *fSimple = new TFile(filenameSimple.c_str(),"RECREATE");
    if(fSimple->IsZombie())
    {
        std::cout << "Cannot create " << filenameSimple << std::endl;
        throw "File is zombie";
    }

    fSimple->cd();

    std::cout<<"Filling Simple Top Tagger Histo ;) "<<std::endl;
    for(int i = 0; i < fakeratePlotsSimple.histos_.size(); i++)
    {
        fakeratePlotsSimple.histos_[i]->Write();
    }

    fSimple->Write();
    fSimple->Close();      

    ///////////////////////////////////
    //Medium Top Tagger
    ///////////////////////////////////  

    const char* rootFileMedium = "joesGroup/TT_TTbarSingleLep_mediumTopTagger.root";    
    std::string filenameMedium = "efficiencyandFakeRatePlots_TT_TTbarSingleLep_mediumTopTagger.root";

    FakeRatePlots fakeratePlotsMedium;    
    fakeratePlotsMedium.makeTH1F("Lep0_", rootFileMedium);
    fakeratePlotsMedium.makeTH1F("Lep1_", rootFileMedium);

    TFile *fMedium = new TFile(filenameMedium.c_str(),"RECREATE");
    if(fMedium->IsZombie())
    {
        std::cout << "Cannot create " << filenameMedium << std::endl;
        throw "File is zombie";
    }

    fMedium->cd();

    std::cout<<"Filling Medium Top Tagger Histo ;) "<<std::endl;
    for(int i = 0; i < fakeratePlotsMedium.histos_.size(); i++)
    {
        fakeratePlotsMedium.histos_[i]->Write();
    }

    fMedium->Write();
    fMedium->Close();

    for(int i = 0; i < fakeratePlotsSimple.histos_.size(); i++)
    {
        makePlots(fakeratePlotsSimple.histoName_[i], fakeratePlotsSimple.histos_[i], fakeratePlotsMedium.histos_[i]);
    }
}
