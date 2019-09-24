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
#include "TSystem.h"

#include <iostream>
#include <string>

#include <cstdlib>
#include <cstdio>
#include <cstring>

class TaggerInfo
{
public:
    std::string inputFile;
    std::string outputFile;
    std::string legName;
    std::string wp;
    
    TaggerInfo(std::string inputFile, std::string outputFile, std::string legName, std::string wp) : inputFile(inputFile), outputFile(outputFile), legName(legName), wp(wp) {} 
    ~TaggerInfo(){}
};

void makePlots(const std::string& outFile, const std::string& name, TH1* simpleHist, TH1* mediumHist,
               const std::string& dataSet, const std::string& simpleLeg = "WP 0.92", const std::string& mediumLeg = "WP 0.99")
{
    //Define canvas and legend
    TCanvas *c = new TCanvas( (dataSet+name).c_str(),(dataSet+name).c_str(),1000,800);  
    TLegend *l = new TLegend(0.68, 0.8, 0.99, 0.9); 
    double topM,     bottomM,      rightM,      leftM;
    /*  */ topM=0.1; bottomM=0.35; rightM=0.05; leftM=0.16;
    gPad->SetAstat(0);
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
    double dw_correction = 1.0;   //40;
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
    l->AddEntry(mediumHist, (mediumLeg).c_str()   , "l");
    l->AddEntry(simpleHist, (simpleLeg).c_str()   , "l");
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

    //std::cout<<"fakeRateEfficiencyPlots/" + dataSet + outFile<<std::endl;
    gSystem -> Exec( ("mkdir -p fakeRateEfficiencyPlots/" + dataSet + outFile).c_str() ) ;    
    c->SaveAs( ( "fakeRateEfficiencyPlots/" + dataSet + name + ".png" ).c_str() );        

    delete l;
    delete pad_up;
    delete pad_dw;
    delete line;
}

void runPlotter(const std::vector<TaggerInfo>& taggerInfo, const std::vector<std::string>& selections, const std::string& dataSet)
{
    char copy[128];
    strcpy(copy, taggerInfo[0].inputFile.c_str() );   
    char* type1;
    type1 = strtok( copy, "-/" );
    char* type2;
    type2 = strtok( nullptr, "-/" );
    char* type3;
    type3 = strtok( nullptr, "-/" );

    // -----------------------------
    // -- Looping over TaggerInfo
    // -----------------------------
    std::vector<Eff_FakeRatePlots*> vecPlots;
    for(const auto& tI : taggerInfo)
    {
        Eff_FakeRatePlots* fakeratePlots = new Eff_FakeRatePlots();
        for(const auto& s : selections)
        {
            fakeratePlots->makeTH1F(s, tI.inputFile.c_str(), type3, tI.wp); 
        }

        TFile *f = new TFile(tI.outputFile.c_str(),"RECREATE");
        if(f->IsZombie())
        {
            std::cout << "Cannot create " << tI.outputFile << std::endl;
            throw "File is zombie";
        }
    
        f->cd();
    
        for(int i = 0; i < fakeratePlots->histos_.size(); i++)
        {
            fakeratePlots->histos_[i]->Write();
        }
    
        f->Write();
        //f->Close();

        vecPlots.push_back(fakeratePlots);
    }
    
    // ------------------
    // -- Making Plots
    // ------------------
    for(int i = 0; i < vecPlots[0]->histos_.size(); i++)
    {
        makePlots(vecPlots[0]->outFile_[i], vecPlots[0]->histoName_[i], vecPlots[0]->histos_[i], vecPlots[1]->histos_[i], dataSet, taggerInfo[0].legName, taggerInfo[1].legName);
    }

    //Cleaning up
    for(auto p : vecPlots)
    {
        delete p;
    }
}

int main()
{
    TH1::AddDirectory(false);
    std::vector<std::string> selections = {"histos/", "Njet7/", "Njet8/", "Njet9/", "Njet10/", "Njet11/", "Njet12/", "Njet13/"};
    
    //// ----------------
    //// -- 2016 TTbar
    //// ----------------
    //std::vector<TaggerInfo> TTbar2016
    //{
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2016_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2016_TT.root", "outputRoot/EfficiencyFakeRatePlots_2016_TTbar_ResolvedTopTagger92.root", "Resolved - 2016_TT - WP 0.92", "0.920"},
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2016_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2016_TT.root", "outputRoot/EfficiencyFakeRatePlots_2016_TTbar_ResolvedTopTagger99.root", "Resolved - 2016_TT - WP 0.99", "0.990"},
    //};
    //runPlotter(TTbar2016, selections, "2016_TT_compareWP-92-99/");

    //// ----------------
    //// -- 2017 TTbar
    //// ----------------
    //std::vector<TaggerInfo> TTbar2017 
    //{
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2017_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2017_TT.root", "outputRoot/EfficiencyFakeRatePlots_2017_TTbar_ResolvedTopTagger92.root", "Resolved - 2017_TT - WP 0.92", "0.920"},
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2017_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2017_TT.root", "outputRoot/EfficiencyFakeRatePlots_2017_TTbar_ResolvedTopTagger99.root", "Resolved - 2017_TT - WP 0.99", "0.990"}
    //};
    //runPlotter(TTbar2017, selections, "2017_TT_compareWP-92-99/");

    //// ----------------
    //// -- 2018 TTbar
    //// ----------------
    //std::vector<TaggerInfo> TTbar2018
    //{
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2018_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2018_TT.root", "outputRoot/EfficiencyFakeRatePlots_2018_TTbar_ResolvedTopTagger92.root", "Resolved - 2018_TT - WP 0.92", "0.920"},
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2018_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2018_TT.root", "outputRoot/EfficiencyFakeRatePlots_2018_TTbar_ResolvedTopTagger99.root", "Resolved - 2018_TT - WP 0.99", "0.990"},
    //};
    //runPlotter(TTbar2018, selections, "2018_TT_compareWP-92-99/");

    //// --------------
    //// -- 2016 QCD 
    //// --------------
    //std::vector<TaggerInfo> QCD2016
    //{
    //  {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2016_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2016_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2016_QCD_ResolvedTopTagger92.root", "Resolved - 2016_QCD - WP 0.92", "0.920"},
    //  {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2016_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2016_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2016_QCD_ResolvedTopTagger99.root", "Resolved - 2016_QCD - WP 0.99", "0.990"},
    //};
    //runPlotter(QCD2016, selections, "2016_QCD_compareWP-92-99/");
     
    //// --------------
    //// -- 2017 QCD
    //// --------------
    //std::vector<TaggerInfo> QCD2017
    //{
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2017_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2017_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2017_QCD_ResolvedTopTagger92.root", "Resolved - 2017_QCD - WP 0.92", "0.920"},
    //    {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2017_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2017_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2017_QCD_ResolvedTopTagger99.root", "Resolved - 2017_QCD - WP 0.99", "0.990"},
    //};
    //runPlotter(QCD2017, selections, "2017_QCD_compareWP-92-99/");

    // --------------
    // -- 2018 QCD
    // --------------
    std::vector<TaggerInfo> QCD2018
    {
        {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2018_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2018_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2018_QCD_ResolvedTopTagger92.root", "Resolved - 2018_QCD - WP 0.92", "0.920"},
        {"/uscms_data/d3/semrat/CMSSW_9_3_3/src/Analyzer/Analyzer/test/condor/hadd_2018_MC_EfficienyCompNjets_ResolvedTTwp0.0_11thSlides_25.09.2019/2018_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2018_QCD_ResolvedTopTagger99.root", "Resolved - 2018_QCD - WP 0.99", "0.990"},
    };
    runPlotter(QCD2018, selections, "2018_QCD_compareWP-92-99/");

}
