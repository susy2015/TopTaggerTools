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

#include <cstdlib>
#include <cstdio>
#include <cstring>

class FakeRatePlots
{

private:  
  
public:
    std::vector<TH1*> histos_;
    std::vector<std::string> histoName_;
    
    TH1F* defineTH1F(const char *file="", std::string branch="")
    {
        //TFile *f = new TFile(file);
        TFile *f = TFile::Open(file);
        TH1F *h = (TH1F*)f->Get(branch.c_str());
        f->Close();
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
  
    void ratioTH1F(std::string fakerateVar, std::string Var, const char* file, std::string fakerate, std::string type, std::string xname, std::string yname)
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
        fakeRate->SetTitle((type + "/" + fakerate).c_str());
        fakeRate->GetXaxis()->SetTitle(xname.c_str());
        fakeRate->GetYaxis()->SetTitle(yname.c_str());
        histos_.emplace_back(fakeRate);
        histoName_.emplace_back(type + "/" + fakerate);
    }

    void makeTH1F(std::string name, const char* rootFile, std::string type)
    {
        ratioTH1F( name + "fakerateMET"     ,name + "MET"               ,rootFile ,name + "fakerateMET"              ,type ,"MET (GeV)"     ,"Fakerate"  );
        ratioTH1F( name + "fakerateNj"      ,name + "nJets"             ,rootFile ,name + "fakerateNjets"            ,type ,"N_{Jets}"      ,"Fakerate"  );
        ratioTH1F( name + "fakerateNb"      ,name + "nBJets"            ,rootFile ,name + "fakerateNbjets"           ,type ,"N_{BJets}"     ,"Fakerate"  );			     
        ratioTH1F( name + "fakerateMET2"    ,name + "MET"               ,rootFile ,name + "fakerateMET2"             ,type ,"MET (GeV)"     ,"Fakerate"  );
        ratioTH1F( name + "fakerateNj2"     ,name + "nJets"             ,rootFile ,name + "fakerateNjets2"           ,type ,"N Jets"        ,"Fakerate"  );
        ratioTH1F( name + "fakerateNb2"     ,name + "nBJets"            ,rootFile ,name + "fakerateNbjets2"          ,type ,"N BJets"       ,"Fakerate"  );    			     
        ratioTH1F( name + "randomTopPt"     ,name + "randomTopCandPt"   ,rootFile ,name + "fakerateRandomTopPt"      ,type ,"P_{T} (GeV)"   ,"Fakerate"  );
        ratioTH1F( name + "randomTopMass"   ,name + "randomTopCandMass" ,rootFile ,name + "fakerateRandomTopMass"    ,type ,"Mass (GeV)"    ,"Fakerate"  );
        ratioTH1F( name + "randomTopEta"    ,name + "randomTopCandEta"  ,rootFile ,name + "fakerateRandomTopEta"     ,type ,"#eta"          ,"Fakerate"  );    
        ratioTH1F( name + "genTopMatchPt"   ,name + "genTopPt"          ,rootFile ,name + "efficiencyGenTopMatchPt"  ,type ,"P_{T} (GeV)"   ,"Efficiency");
        ratioTH1F( name + "genTopMatchMass" ,name + "genTopMass"        ,rootFile ,name + "efficiencyGenTopMatchMass",type ,"Mass (GeV)"    ,"Efficiency");	
        ratioTH1F( name + "genTopMatchEta"  ,name + "genTopEta"         ,rootFile ,name + "efficiencyGenTopMatchEta" ,type ,"#eta"          ,"Efficiency");
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
    
    mediumHist->SetMaximum( 1.2*( mediumHist->GetMaximum() ) );
    mediumHist->SetMinimum(0);
    mediumHist->SetTitleSize(0.002);
    mediumHist->SetTitleSize(0.05,"X");
    mediumHist->SetTitleSize(0.05,"Y");
    mediumHist->SetTitleOffset(1.2,"X");
    mediumHist->SetTitleOffset(1.5,"Y");
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

void runPlotter(const char* rootFileSimple, std::string filenameSimple, const char* rootFileMedium, std::string filenameMedium)
{

    char copy[128];
    strcpy(copy, rootFileSimple);    
    char* type1;
    type1 = strtok( copy, "-/" );
    char* type2;
    type2 = strtok( nullptr, "-/" );
    char* type3;
    type3 = strtok( nullptr, "-/" );
    
    ///////////////////////////////
    //      Simple TopTagger
    ///////////////////////////////
    FakeRatePlots fakeratePlotsSimple;
    fakeratePlotsSimple.makeTH1F("Lep0/", rootFileSimple, type3);
    fakeratePlotsSimple.makeTH1F("Lep1/", rootFileSimple, type3);

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
    
    ///////////////////////////////
    //      Medium TopTagger
    ///////////////////////////////    
    FakeRatePlots fakeratePlotsMedium;
    fakeratePlotsMedium.makeTH1F("Lep0/", rootFileMedium, type3);
    fakeratePlotsMedium.makeTH1F("Lep1/", rootFileMedium, type3);
    
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
    
    ////////////////////////////////
    //       Making Plots
    ////////////////////////////////
    for(int i = 0; i < fakeratePlotsSimple.histos_.size(); i++)
    {
        makePlots(fakeratePlotsSimple.histoName_[i], fakeratePlotsSimple.histos_[i], fakeratePlotsMedium.histos_[i]);
    }
}

int main()
{
    TH1::AddDirectory(false);

    runPlotter("joesGroup/SimpleHaddFiles/TT_TTbar-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbar_simpleTopTagger.root",
               "joesGroup/MediumHaddFiles/TT_TTbar-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbar_mediumTopTagger.root"
              );
    runPlotter("joesGroup/SimpleHaddFiles/TT_QCD-2018-3-9.root", "efficiencyandFakeRatePlots_TT_QCD_simpleTopTagger.root",
               "joesGroup/MediumHaddFiles/TT_QCD-2018-3-9.root", "efficiencyandFakeRatePlots_TT_QCD_mediumTopTagger.root"
              );
    runPlotter("joesGroup/SimpleHaddFiles/TT_TTbarSingleLep-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root",
               "joesGroup/MediumHaddFiles/TT_TTbarSingleLep-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbarSingleLep_mediumTopTagger.root"
              );
}
