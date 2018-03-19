//#include "../include/Eff_FakeRatePlots.h"
#include "Eff_FakeRatePlots.h"

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

void makePlots(std::string name, TH1* simpleHist, TH1* mediumHist, std::string dataSet)
{
    //Define canvas and legend
    TCanvas *c = new TCanvas( (dataSet+name).c_str(),(dataSet+name).c_str(),1000,800);  
    TLegend *l = new TLegend(0.68, 0.8, 0.99, 0.9);
    double topM,     bottomM,      rightM,      leftM;
    /*  */ topM=0.1; bottomM=0.35; rightM=0.05; leftM=0.16;
    gPad->SetTopMargin(   0.0);
    gPad->SetBottomMargin(0.0);
    gPad->SetRightMargin( 0.0);
    gPad->SetLeftMargin(  0.0);
    //gPad->SetTopMargin(topM);
    //gPad->SetBottomMargin(bottomM);
    //gPad->SetRightMargin(rightM);
    //gPad->SetLeftMargin(leftM);

    //Define the top and bottom TPad
    double up_height     = 0.75;  // please tune so that the upper figures size will meet your requirement
    double dw_correction = 1.0;//40;
    double dw_height     = (1.0 - up_height)*dw_correction;
    TPad *pad_up = new TPad("pad_up1","pad_up1",0.0, 1.0 - up_height, 1.0,       1.0);
    TPad *pad_dw = new TPad("pad_dw1","pad_dw1",0.0, 0.0,             1.0, dw_height);
    pad_up->Draw();
    pad_dw->Draw();
  
    pad_up->SetFrameFillColor(0);
    pad_up->SetFillColor(0);
    pad_up->SetTopMargin(topM);
    pad_up->SetBottomMargin(0.02);
    pad_up->SetLeftMargin(leftM);
    pad_up->SetRightMargin(rightM);

    pad_dw->SetFrameFillColor(0);
    pad_dw->SetFillColor(0);
    pad_dw->SetTopMargin(0.04);
    pad_dw->SetBottomMargin(bottomM);
    pad_dw->SetLeftMargin(leftM);
    pad_dw->SetRightMargin(rightM);

    //Top TPad
    pad_up->cd();
    
    double max = std::max( mediumHist->GetMaximum(), simpleHist->GetMaximum() );
    mediumHist->SetMaximum( 1.2*max );
    mediumHist->SetMinimum(0);
    mediumHist->SetTitleSize(0.002);
    mediumHist->SetTitleSize(0.05,"X");
    mediumHist->SetTitleSize(0.05,"Y");
    mediumHist->SetTitleOffset(1.2,"X");
    mediumHist->SetTitleOffset(1.5,"Y");
    mediumHist->SetLabelSize(0.0,"X");
    //mediumHist->SetLabelSize(0.05,"X");
    mediumHist->SetLabelSize(0.05,"Y");
    //mediumHist->SetNdivisions(0,"X");
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

    //Bottom TPad
    pad_dw->cd();

    Eff_FakeRatePlots dummy;
    TH1F *ratio = (TH1F*)mediumHist->Clone( (name+dataSet+"Ratio").c_str() );;
    for(int i=0; i < ratio->GetSize(); i++)
    {
        double errorNum = dummy.errorRatio(i,mediumHist,simpleHist);
        double ratioNum = mediumHist->GetBinContent(i)/simpleHist->GetBinContent(i);      
        if(isinf(ratioNum)==1 || isnan(ratioNum)==1) ratioNum = 0;
        ratio->SetBinContent(i,ratioNum);
        ratio->SetBinError(i,errorNum);
    }

    double maxDw = ratio->GetMaximum();
    double minDw = ratio->GetMinimum();

    ratio->SetTitle("");
    ratio->SetStats(false);
    ratio->SetMaximum( 1.5 );
    ratio->SetMinimum( 0.5 );
    //ratio->SetMaximum( std::max(1.5*maxDw, 1.3) );
    //ratio->SetMinimum( std::min(0.5*maxDw, 0.7) );
    ratio->SetTitleOffset(1.0,"X");
    ratio->SetTitleOffset(0.5,"Y");
    ratio->SetTitleSize(0.15,"X");
    ratio->SetTitleSize(0.15,"Y");
    ratio->SetLabelSize(0.15,"X");
    ratio->SetLabelSize(0.15,"Y");
    ratio->SetNdivisions(10,"X");
    ratio->SetNdivisions( 5,"Y");
    ratio->GetYaxis()->SetTitle("Ratio");
    ratio->SetLineColor(kBlack);
    ratio->Draw("E1");

    TF1 *line = new TF1( (name+dataSet+"Line").c_str(),"1",-2000,2000);
    line->SetLineColor(kRed);
    line->Draw("same");
    
    c->SaveAs( ("plots/" + dataSet + name + ".png").c_str() );        
}

void runPlotter(const char* rootFileSimple, std::string filenameSimple, const char* rootFileMedium, std::string filenameMedium, std::string dataSet)
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
    Eff_FakeRatePlots fakeratePlotsSimple;
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
    Eff_FakeRatePlots fakeratePlotsMedium;
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
        makePlots(fakeratePlotsSimple.histoName_[i], fakeratePlotsSimple.histos_[i], fakeratePlotsMedium.histos_[i], dataSet);
    }
}

int main()
{
    TH1::AddDirectory(false);

    ////Orginal cuts for Simple and Medium Top Tagger
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbar-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbar_simpleTopTagger.root",
    //           "joesGroup/MediumHaddFiles/TT_TTbar-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbar_mediumTopTagger.root",
    //           "originalDiscCuts/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_QCD-2018-3-9.root", "efficiencyandFakeRatePlots_TT_QCD_simpleTopTagger.root",
    //           "joesGroup/MediumHaddFiles/TT_QCD-2018-3-9.root", "efficiencyandFakeRatePlots_TT_QCD_mediumTopTagger.root",
    //           "originalDiscCuts/"               
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbarSingleLep-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root",
    //           "joesGroup/MediumHaddFiles/TT_TTbarSingleLep-2018-3-9.root", "efficiencyandFakeRatePlots_TT_TTbarSingleLep_mediumTopTagger.root",
    //           "originalDiscCuts/"               
    //          );
    //
    //
    ////95max_00017pt85 Medium Top Tagger cuts
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbar-2018-3-9.root"                 ,"efficiencyandFakeRatePlots_TT_TTbar_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_00017pt85/TT_TTbar-2018-3-13.root","efficiencyandFakeRatePlots_TT_TTbar_95max_00017pt85_medium.root",
    //           "95max_00017pt85_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_QCD-2018-3-9.root"                 ,"efficiencyandFakeRatePlots_TT_QCD_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_00017pt85/TT_QCD-2018-3-13.root","efficiencyandFakeRatePlots_TT_QCD_95max_00017pt85_medium.root",
    //           "95max_00017pt85_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbarSingleLep-2018-3-9.root"                 ,"efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_00017pt85/TT_TTbarSingleLep-2018-3-13.root","efficiencyandFakeRatePlots_TT_TTbarSingleLep_95max_00017pt85_medium.root",
    //           "95max_00017pt85_medium/"
    //          );
    //
    //
    ////95max_0004375pt775 Medium Top Tagger cuts
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbar-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbar_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_0004375pt775/TT_TTbar-2018-3-12.root","efficiencyandFakeRatePlots_TT_TTbar_95max_0004375pt775_medium.root",
    //           "95max_0004375pt775_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_QCD-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_QCD_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_0004375pt775/TT_QCD-2018-3-12.root","efficiencyandFakeRatePlots_TT_QCD_95max_0004375pt775_medium.root",
    //           "95max_0004375pt775_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbarSingleLep-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_0004375pt775/TT_TTbarSingleLep-2018-3-12.root","efficiencyandFakeRatePlots_TT_TTbarSingleLep_95max_0004375pt775_medium.root",
    //           "95max_0004375pt775_medium/"
    //          );
    //
    ////95max_0005pt7 Medium Top Tagger cuts
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbar-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbar_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_0005pt7/TT_TTbar-2018-3-12.root","efficiencyandFakeRatePlots_TT_TTbar_95max_0005pt7_medium.root",
    //           "95max_0005pt7_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_QCD-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_QCD_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_0005pt7/TT_QCD-2018-3-12.root","efficiencyandFakeRatePlots_TT_QCD_95max_0005pt7_medium.root",
    //           "95max_0005pt7_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbarSingleLep-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_0005pt7/TT_TTbarSingleLep-2018-3-12.root","efficiencyandFakeRatePlots_TT_TTbarSingleLep_95max_0005pt7_medium.root",
    //           "95max_0005pt7_medium/"
    //          );
    //
    ////95max_00175pt3125 Medium Top Tagger cuts
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbar-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbar_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_00175pt3125/TT_TTbar-2018-3-13.root","efficiencyandFakeRatePlots_TT_TTbar_95max_00175pt3125_medium.root",
    //           "95max_00175pt3125_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_QCD-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_QCD_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_00175pt3125/TT_QCD-2018-3-13.root","efficiencyandFakeRatePlots_TT_QCD_95max_00175pt3125_medium.root",
    //           "95max_00175pt3125_medium/"
    //          );
    //runPlotter("joesGroup/SimpleHaddFiles/TT_TTbarSingleLep-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root",
    //           "joesGroup/mediumTopTagger_95max_00175pt3125/TT_TTbarSingleLep-2018-3-13.root","efficiencyandFakeRatePlots_TT_TTbarSingleLep_95max_00175pt3125_medium.root",
    //           "95max_00175pt3125_medium/"
    //          );    

    //93max_000575pt7 Medium Top Tagger cuts
    runPlotter("joesGroup/SimpleHaddFiles/TT_TTbar-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbar_simpleTopTagger.root",
               "joesGroup/mediumTopTagger_93max_000575pt7/TT_TTbar-2018-3-13.root","efficiencyandFakeRatePlots_TT_TTbar_93max_000575pt7_medium.root",
               "93max_000575pt7_medium/"
              );
    runPlotter("joesGroup/SimpleHaddFiles/TT_QCD-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_QCD_simpleTopTagger.root",
               "joesGroup/mediumTopTagger_93max_000575pt7/TT_QCD-2018-3-13.root","efficiencyandFakeRatePlots_TT_QCD_93max_000575pt7_medium.root",
               "93max_000575pt7_medium/"
              );
    runPlotter("joesGroup/SimpleHaddFiles/TT_TTbarSingleLep-2018-3-9.root"                    ,"efficiencyandFakeRatePlots_TT_TTbarSingleLep_simpleTopTagger.root",
               "joesGroup/mediumTopTagger_93max_000575pt7/TT_TTbarSingleLep-2018-3-13.root","efficiencyandFakeRatePlots_TT_TTbarSingleLep_93max_000575pt7_medium.root",
               "93max_000575pt7_medium/"
              );    
}
