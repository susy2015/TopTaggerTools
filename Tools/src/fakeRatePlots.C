#include "Eff_FakeRatePlots.h"

#include "TColor.h"
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

#include <getopt.h>
#include <unistd.h>

class TaggerInfo
{
public:
    std::string inputFile;
    std::string outputFile;
    std::string legName;
    std::string wp;
    std::string dataset;
 
    TaggerInfo(std::string inputFile, std::string outputFile, std::string legName, std::string wp, std::string dataset) : inputFile(inputFile), outputFile(outputFile), legName(legName), wp(wp), dataset(dataset){} 
    ~TaggerInfo(){}
};

void makePlots(const std::string& outFile, const std::string& name, TH1* simpleHist, TH1* mediumHist,
               const std::string& dataSet, const std::string& simpleLeg = "WP 0.92", const std::string& mediumLeg = "WP 0.95", 
               const std::string& year = "2016")
{
    //Define canvas and legend
    TCanvas *c = new TCanvas( (dataSet+name).c_str(),(dataSet+name).c_str(),800,800);  
    TLegend *l = new TLegend(0.30, 0.77, 0.99, 0.87); // 0.68, 0.8, 0.99, 0.9 
    double topM,     bottomM,      rightM,      leftM;
           topM=0.1; bottomM=0.35; rightM=0.05; leftM=0.14;
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
    double dw_correction = 1.0;   // 40;
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
    pad_up->SetGrid();
 
    double max = std::max( mediumHist->GetMaximum(), simpleHist->GetMaximum() );
    mediumHist->SetMaximum( 1.2*max );
    mediumHist->SetMinimum(0);
    mediumHist->SetTitle("");
    mediumHist->SetTitleSize(0.002);
    mediumHist->SetTitleSize(0.05,"X");
    mediumHist->SetTitleSize(0.05,"Y");
    mediumHist->SetTitleOffset(1.2,"X");
    mediumHist->SetTitleOffset(0.5 * (0.75/0.25),"Y");
    mediumHist->SetLabelSize(0.0,"X");
    mediumHist->SetLabelSize(0.05,"Y");
    mediumHist->SetStats(false);   
    mediumHist->SetLineColor(kCyan+2);
    mediumHist->SetLineWidth(2);
    mediumHist->Draw("hist E");

    simpleHist->SetLineColor(kMagenta+2);
    simpleHist->SetLineWidth(2);
    simpleHist->Draw("hist E same");

    TLatex mark;
    mark.SetNDC(true);
    mark.SetTextAlign(11);
    mark.SetTextSize(0.06);
    mark.SetTextFont(61);
    mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS");
    mark.SetTextSize(0.040);
    mark.SetTextFont(52);
    mark.DrawLatex(gPad->GetLeftMargin() + 0.1, 1 - (gPad->GetTopMargin() - 0.017), "Work in Progress");
    mark.SetTextSize(0.040);
    mark.SetTextFont(42);
    mark.SetTextAlign(31);
    mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), (year + " (13 TeV)").c_str());

    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(0.03);    
    l->AddEntry(mediumHist, (mediumLeg).c_str(), "l");
    l->AddEntry(simpleHist, (simpleLeg).c_str(), "l");
    l->Draw();

    //Bottom TPad
    pad_dw->cd();
    pad_dw->SetGrid();

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
    ratio->SetMaximum( 1.75 ); 
    ratio->SetMinimum( 0.25 ); 
    ratio->SetTitleOffset(1.0,"X");
    ratio->SetTitleOffset(0.3,"Y");
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
    line->SetLineColor(kBlue);
    line->SetLineWidth(1);
    line->Draw("same");

    gSystem -> Exec( ("mkdir -p fakeRateEfficiencyPlots/" + dataSet + outFile).c_str() ) ;    
    c->SaveAs( ( "fakeRateEfficiencyPlots/" + dataSet + name + ".pdf" ).c_str() );        

    delete l;
    delete pad_up;
    delete pad_dw;
    delete line;
}

void makePlotsThesis(const std::string& outFile, const std::string& name, TH1* hist1, TH1* hist2, TH1* hist3, TH1* hist4, TH1* hist5,
               const std::string& dataSet, const std::string& leg1 = "WP 0.92", const std::string& leg2 = "WP 0.95", 
               const std::string& leg3 = "WP 0.96", const std::string& leg4 = "WP 0.97", 
               const std::string& leg5 = "WP 0.98", const std::string& year = "2016")
{
    //Define canvas and legend
    TCanvas *c = new TCanvas( (dataSet+name).c_str(),(dataSet+name).c_str(),800,800);  
    TLegend *l = new TLegend(0.30, 0.67, 0.99, 0.87); // 0.68, 0.8, 0.99, 0.9 
    double topM,     bottomM,      rightM,      leftM;
           topM=0.1; bottomM=0.35; rightM=0.05; leftM=0.11;
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
    double up_height     = 1.0;  // please tune so that the upper figures size will meet your requirement
    double dw_correction = 1.0;   // 40;
    double dw_height     = (1.0 - up_height)*dw_correction;
    TPad *pad_up = new TPad("pad_up1","pad_up1",0.0, 1.0 - up_height, 1.0,       1.0);
    TPad *pad_dw = new TPad("pad_dw1","pad_dw1",0.0, 0.0,             1.0, dw_height);
    pad_up->Draw();
    pad_dw->Draw();
    pad_up->SetFrameFillColor(0);
    pad_up->SetFillColor(0);
    pad_up->SetTopMargin(topM);
    pad_up->SetBottomMargin(0.12);
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
    pad_up->SetGrid();
    double max = std::max( hist1->GetMaximum(), hist2->GetMaximum() );
    hist1->SetMaximum( 1.2*max );
    hist1->SetMinimum(0);
    hist1->SetTitle("");
    hist1->SetTitleSize(0.002);
    hist1->SetTitleSize(0.05,"X");
    hist1->SetTitleSize(0.05,"Y");
    hist1->SetTitleOffset(1.0,"X");
    hist1->SetTitleOffset(0.35 * (0.75/0.25),"Y");
    hist1->SetLabelSize(0.035,"X");
    hist1->SetLabelSize(0.035,"Y");
    hist1->SetStats(false);   
    hist1->SetLineColor(TColor::GetColor("#cd5c5c"));
    hist2->SetLineColor(TColor::GetColor("#b2e2e2"));
    hist3->SetLineColor(TColor::GetColor("#66c2a4"));
    hist4->SetLineColor(TColor::GetColor("#2ca25f"));
    hist5->SetLineColor(TColor::GetColor("#00441b"));
    //hist5->SetLineColor(TColor::GetColor("#006d2c"));

    hist1->SetLineWidth(2);
    hist1->Draw("hist E");

    hist2->SetLineWidth(2);
    hist2->Draw("hist E same");

    hist3->SetLineWidth(2);
    hist3->Draw("hist E same");

    hist4->SetLineWidth(2);
    hist4->Draw("hist E same");

    hist5->SetLineWidth(2);
    hist5->Draw("hist E same");

    hist1->Draw("hist E same");

    TLatex mark;
    mark.SetNDC(true);
    mark.SetTextAlign(11);
    mark.SetTextSize(0.06);
    mark.SetTextFont(61);
    mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS");
    mark.SetTextSize(0.040);
    mark.SetTextFont(52);
    //mark.DrawLatex(gPad->GetLeftMargin() + 0.1, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");
    mark.DrawLatex(gPad->GetLeftMargin() + 0.1, 1 - (gPad->GetTopMargin() - 0.017), "   Work in Progress");
    mark.SetTextSize(0.040);
    mark.SetTextFont(42);
    mark.SetTextAlign(31);
    mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), (year + " (13 TeV)").c_str());

    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextSize(0.03);    
    l->AddEntry(hist1, (leg1).c_str(), "l");
    l->AddEntry(hist2, (leg2).c_str(), "l");
    l->AddEntry(hist3, (leg3).c_str(), "l");
    l->AddEntry(hist4, (leg4).c_str(), "l");
    l->AddEntry(hist5, (leg5).c_str(), "l");
    l->Draw();

    //Bottom TPad
    //pad_dw->cd();
    //pad_dw->SetGrid();

    //Eff_FakeRatePlots dummy;
    //TH1F *ratio = (TH1F*)mediumHist->Clone( (name+dataSet+"Ratio").c_str() );;
    //for(int i=0; i < ratio->GetSize(); i++)
    //{
    //    double errorNum = dummy.errorRatio(i,mediumHist,simpleHist); 
    //    double ratioNum = mediumHist->GetBinContent(i)/simpleHist->GetBinContent(i);      
    //    if(isinf(ratioNum)==1 || isnan(ratioNum)==1) ratioNum = 0;
    //    ratio->SetBinContent(i,ratioNum);
    //    ratio->SetBinError(i,errorNum);
    //}
    //double maxDw = ratio->GetMaximum();
    //double minDw = ratio->GetMinimum();

    //ratio->SetTitle("");
    //ratio->SetStats(false);
    //ratio->SetMaximum( 1.75 ); 
    //ratio->SetMinimum( 0.25 ); 
    //ratio->SetTitleOffset(1.0,"X");
    //ratio->SetTitleOffset(0.3,"Y");
    //ratio->SetTitleSize(0.15,"X");
    //ratio->SetTitleSize(0.15,"Y");
    //ratio->SetLabelSize(0.15,"X");
    //ratio->SetLabelSize(0.15,"Y");
    //ratio->SetNdivisions(10,"X");
    //ratio->SetNdivisions( 5,"Y");
    //ratio->GetYaxis()->SetTitle("Ratio");
    //ratio->SetLineColor(kBlack);
    //ratio->Draw("E1");

    //TF1 *line = new TF1( (name+dataSet+"Line").c_str(),"1",-2000,2000);
    //line->SetLineColor(kBlue);
    //line->SetLineWidth(1);
    //line->Draw("same");

    gSystem -> Exec( ("mkdir -p fakeRateEfficiencyPlotsThesis/" + dataSet + outFile).c_str() ) ;    
    c->SaveAs( ( "fakeRateEfficiencyPlotsThesis/" + dataSet + name + ".pdf" ).c_str() );        

    delete l;
    delete pad_up;
    delete pad_dw;
    //delete line;
}


void runPlotter(const std::vector<TaggerInfo>& taggerInfo, const std::vector<std::string>& selections, 
                const std::string& dataSet, const std::string& year)
{
    char copy[256]; // 128
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
            fakeratePlots->makeTH1F(s, tI.inputFile.c_str(), type3, tI.wp, tI.dataset);
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
        if (vecPlots.size() < 3) {
            makePlots(vecPlots[0]->outFile_[i], vecPlots[0]->histoName_[i], vecPlots[0]->histos_[i], vecPlots[1]->histos_[i], dataSet, taggerInfo[0].legName, taggerInfo[1].legName, year);
        } else
        {
            makePlotsThesis(vecPlots[0]->outFile_[i], vecPlots[0]->histoName_[i], vecPlots[0]->histos_[i], vecPlots[1]->histos_[i],vecPlots[2]->histos_[i], vecPlots[3]->histos_[i],vecPlots[4]->histos_[i], dataSet, taggerInfo[0].legName, taggerInfo[1].legName, taggerInfo[2].legName, taggerInfo[3].legName, taggerInfo[4].legName, year);
        }

    }

    //Cleaning up
    for(auto p : vecPlots)
    {
        delete p;
    }
}

int main(int argc, char *argv[])
{
    int opt, option_index = 0;
    std::string year = "", dataset = "", wp1 = "", wp2 = "";
    static struct option long_options[] = {
        {"year",    required_argument, 0, 'y'},
        {"dataset", required_argument, 0, 'd'},
        {"wp1",     required_argument, 0, 'w'},
        {"wp2",     required_argument, 0, 'x'},
    };

    while((opt = getopt_long(argc, argv, "y:d:w:x:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'y': year    = optarg; break;
            case 'd': dataset = optarg; break;
            case 'w': wp1     = optarg; break;
            case 'x': wp2     = optarg; break;
        }
    }

    TH1::AddDirectory(false);
    std::vector<std::string> selections = {"histos"};
    //std::vector<std::string> selections = {"histos", "Njet7", "Njet8", "Njet9", "Njet10", "Njet11", "Njet12", "Njet12inc"};

    // -------------------
    // -- 2016preVFP TTbar
    // -------------------
    std::vector<TaggerInfo> TTbar2016preVFP
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_TT.root", "outputRoot/EfficiencyFakeRatePlots_2016preVFP_TTbar_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_TT.root", "outputRoot/EfficiencyFakeRatePlots_2016preVFP_TTbar_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset},
    };

    std::vector<TaggerInfo> TTbar2016preVFPThesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_TT.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_TT.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_TT.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_TT.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_TT.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };

    // --------------------
    // -- 2016postVFP TTbar
    // --------------------
    std::vector<TaggerInfo> TTbar2016postVFP
    {                                 
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_TT.root", "outputRoot/EfficiencyFakeRatePlots_2016postVFP_TTbar_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_TT.root", "outputRoot/EfficiencyFakeRatePlots_2016postVFP_TTbar_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset},
    };

    std::vector<TaggerInfo> TTbar2016postVFPThesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_TT.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_TT.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_TT.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_TT.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_TT.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };

    // ----------------
    // -- 2017 TTbar
    // ----------------
    std::vector<TaggerInfo> TTbar2017 
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_TT.root", "outputRoot/EfficiencyFakeRatePlots_2017_TTbar_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset}, 
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_TT.root", "outputRoot/EfficiencyFakeRatePlots_2017_TTbar_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset},
    };

    std::vector<TaggerInfo> TTbar2017Thesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_TT.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_TT.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_TT.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_TT.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_TT.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };

    // -------------------
    // -- 2018 TTbar
    // -------------------
    std::vector<TaggerInfo> TTbar2018
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_TT.root", "outputRoot/EfficiencyFakeRatePlots_2018_TT_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_TT.root", "outputRoot/EfficiencyFakeRatePlots_2018_TT_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset},        
    };

    std::vector<TaggerInfo> TTbar2018Thesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_TT.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_TT.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_TT.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_TT.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_TT.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };

    // -----------------
    // -- 2016preVFP QCD 
    // -----------------
    std::vector<TaggerInfo> QCD2016preVFP
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2016preVFP_QCD_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2016preVFP_QCD_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset},
    };

    std::vector<TaggerInfo> QCD2016preVFPThesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_QCD.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_QCD.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_QCD.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_QCD.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016preVFP_QCD.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };

    // ------------------
    // -- 2016postVFP QCD 
    // ------------------
    std::vector<TaggerInfo> QCD2016postVFP
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2016postVFP_QCD_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2016postVFP_QCD_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset},
    };

    std::vector<TaggerInfo> QCD2016postVFPThesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_QCD.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_QCD.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_QCD.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_QCD.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2016postVFP_QCD.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };


    // --------------
    // -- 2017 QCD
    // --------------
    std::vector<TaggerInfo> QCD2017
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2017_QCD_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2017_QCD_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset}, 
    };

    std::vector<TaggerInfo> QCD2017Thesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_QCD.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_QCD.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_QCD.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_QCD.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2017_QCD.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };

    // -----------------
    // -- 2018 QCD
    // -----------------
    std::vector<TaggerInfo> QCD2018
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2018_QCD_ResolvedTopTagger.root", "WP_" + wp1, wp1, dataset}, 
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_QCD.root", "outputRoot/EfficiencyFakeRatePlots_2018_QCD_ResolvedTopTagger.root", "WP_" + wp2, wp2, dataset},
    };

    std::vector<TaggerInfo> QCD2018Thesis
    {
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_QCD.root", "outputRoot/junk.root", "WP_0.92", "0.920", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_QCD.root", "outputRoot/junk.root", "WP_0.95", "0.950", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_QCD.root", "outputRoot/junk.root", "WP_0.96", "0.960", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_QCD.root", "outputRoot/junk.root", "WP_0.97", "0.970", dataset},
        {"/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_Run2UL_ResolvedTopTagger_fakeRateEfficiency_withDeepFlavor/2018_QCD.root", "outputRoot/junk.root", "WP_0.98", "0.980", dataset},

    };

    // 2016preVFP TT & QCD
    if (year == "2016preVFP" && dataset == "TT")
        { runPlotter(TTbar2016preVFP, selections, "2016preVFP_TT_compareWPs_" + wp1 + "_" + wp2 + "/", year); 
          runPlotter(TTbar2016preVFPThesis, selections, "2016preVFP_TT_compareWPs_Thesis/", year); 
        }

    else if (year == "2016preVFP" && dataset == "QCD")
        { runPlotter(QCD2016preVFP, selections, "2016preVFP_QCD_compareWPs_" + wp1 + "_" + wp2 + "/", year); 
          runPlotter(QCD2016preVFPThesis, selections, "2016preVFP_QCD_compareWPs_Thesis/", year); 
        }

    // 2016postVFP TT & QCD
    if (year == "2016postVFP" && dataset == "TT")
        { runPlotter(TTbar2016postVFP, selections, "2016postVFP_TT_compareWPs_" + wp1 + "_" + wp2 + "/", year); 
          runPlotter(TTbar2016postVFPThesis, selections, "2016postVFP_TT_compareWPs_Thesis/", year); 
        }

    else if (year == "2016postVFP" && dataset == "QCD")
        { runPlotter(QCD2016postVFP, selections, "2016postVFP_QCD_compareWPs_" + wp1 + "_" + wp2 + "/", year);
          runPlotter(QCD2016postVFPThesis, selections, "2016postVFP_QCD_compareWPs_Thesis/", year);
        }

    // 2017 TT & QCD
    else if (year == "2017" && dataset == "TT")
        { runPlotter(TTbar2017, selections, "2017_TT_compareWPs_" + wp1 + "_" + wp2 + "/", year);
          runPlotter(TTbar2017Thesis, selections, "2017_TT_compareWPs_Thesis/", year);
        }

    else if (year == "2017" && dataset == "QCD")
        { runPlotter(QCD2017, selections, "2017_QCD_compareWPs_" + wp1 + "_" + wp2 + "/", year);
          runPlotter(QCD2017Thesis, selections, "2017_QCD_compareWPs_Thesis/", year);
        }

    // 2018 TT & QCD
    else if (year == "2018" && dataset == "TT")
        { runPlotter(TTbar2018, selections, "2018_TT_compareWPs_" + wp1 + "_" + wp2 + "/", year); 
          runPlotter(TTbar2018Thesis, selections, "2018_TT_compareWPs_Thesis/", year); 
        }

    else if (year == "2018" && dataset == "QCD")
        { runPlotter(QCD2018, selections, "2018_QCD_compareWPs_" + wp1 + "_" + wp2 + "/", year);
          runPlotter(QCD2018Thesis, selections, "2018_QCD_compareWPs_Thesis/", year);
        }  

}
