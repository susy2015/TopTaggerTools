##############################################
####  ROC Plots for Resolved Top Tagger:  ####
# (1) Normilize the discriminator histograms #
# (2) Take the integral of them              #
# TPR : True Positive Rate (TTbar)           #
# FPR : False Positive Rate (QCD)            #
##############################################

import ROOT
import math
import sys
import numpy as np
import array as arr
from ROOT import TFile, gROOT, gStyle, TLatex

debug = False 

def print_db(input):
    if (debug):
        print input

def main():
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)   
 
    # ---------------------------------
    # root path & years & histograms
    # --------------------------------- 
    basepath = "/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_YEAR_ResolvedTopTagger_fakeRateEfficiency_withDeepCSV_17.11.2022/"
    
    years = [
        "2016preVFP" ,
        "2016postVFP",
        "2017" ,
        "2018" ,
    ]
    
    histnames = {
        "TT"  : "topDiscGenMatchWP_0.000histos" ,   
        "QCD" : "topDiscNotGenMatchWP_0.000histos" , 
    }
    
    colors = {
        "TT"  : 6 ,
        "QCD" : 7 ,
    }
    
    WPs = [
        0.92 ,
        0.95 ,
        0.96 ,
        0.97 ,
        0.98 ,
        #0.99 ,
    ]
    
    WP_colors = {
        0.92 : ROOT.TColor.GetColor("#cd5c5c"), #46 ,
        0.95 : ROOT.TColor.GetColor("#b2e2e2"), #40 ,
        0.96 : ROOT.TColor.GetColor("#66c2a4"), #38 ,
        0.97 : ROOT.TColor.GetColor("#2ca25f"), #43 ,
        0.98 : ROOT.TColor.GetColor("#00441b"), #30 ,
        #0.99 : 41 ,
    }
    
    label_size     = 0.033
    label_offset_x = 0.03 # 0.02
    label_offset_y = 0.05 # 0.04

    # ----------------------
    # loop over the years 
    # ----------------------
    for year in years:
        print_db("Processing year " + year)
        hists = {}        

        # ------------------------------------------------
        # Normalize & save the discriminator histograms
        # ------------------------------------------------     
        for histname in histnames:
            print_db("Processing type " + histname)
            filename = basepath.replace("YEAR", year) + year + "_" + histname + ".root"
            print_db("Opening file " + filename)
            f = ROOT.TFile.Open(filename, "READ")
            print_db("Getting histogram " + filename + ":histos/" + histnames[histname])
            hists[histname] = f.Get("histos/" + histnames[histname])
            hists[histname].SetDirectory(gROOT)
            hists[histname].SetTitle("")
            print_db("Normalizing histogram " + filename +":histos/" + histnames[histname])
            scale = 1/(hists[histname].Integral())
            hists[histname].Scale(scale)
            f.Close()

        # -------------------- 
        # Get the ROC plots
        # --------------------
        arrays = {
            "TT" : [] ,
            "QCD" : [] ,
        }
        sizes = {}
        for hist in hists:
            for i in range(0, hists[hist].GetNbinsX()):
                arrays[hist].append(hists[hist].Integral(i, hists[hist].GetNbinsX()))
            sizes[hist] = np.array(arrays[hist]).size

        size = 0
        if (sizes["TT"] == sizes["QCD"]):
            size = sizes["TT"]
        else:
            print "ERROR: Size of TT and QCD histograms do not match for year " + year + "!!!"
        
        roc = ROOT.TGraph(size, np.array(arrays["QCD"]).astype(np.double), np.array(arrays["TT"]).astype(np.double))
        print_db(roc)

        # cnavas & legend
        ROC = ROOT.TCanvas("ROC", "ROC", 0, 0, 800, 800)
        ROC.SetGrid()

        roc.SetTitle("")
        #roc.SetTitle(year + " ROC plot")
        roc.GetXaxis().SetTitle("FPR (QCD)")
        roc.GetYaxis().SetTitle("TPR (t#bar{t})")
        roc.Draw("AL")

        # --------------
        # Get the WPs 
        # --------------
        points = {}
        labels = {}
        for WP in WPs:
            bin_num_TT  = hists["TT"].GetXaxis().FindBin(WP)
            bin_num_QCD = hists["QCD"].GetXaxis().FindBin(WP)
            bin_num=0
            if (bin_num_TT == bin_num_QCD):
                bin_num = bin_num_TT
            else:
                print "ERROR: bin_num_TT =/= bin_num_QCD"
            points[WP] = ROOT.TMarker(arrays["QCD"][bin_num],arrays["TT"][bin_num],22)
            #print efficiency and fake rate values for AN
            print "Year -------- : ", (year)
            print "WP        : ", (WP)
            print "Efficiency: ", (arrays["TT"][bin_num])
            print "Fake Rate : ", (arrays["QCD"][bin_num])
            points[WP].SetMarkerSize(3) # 2
            points[WP].SetMarkerColor(WP_colors[WP])
            points[WP].Draw("SAME")
            labels[WP] = ROOT.TLatex(arrays["QCD"][bin_num]+label_offset_x,arrays["TT"][bin_num]-label_offset_y,str(WP))
            labels[WP].SetTextColor(WP_colors[WP])
            labels[WP].SetTextSize(label_size)
            labels[WP].Draw("SAME")

        mark = ROOT.TLatex()
        mark.SetNDC(True)
        mark.SetTextAlign(11)
        mark.SetTextSize(0.06)
        mark.SetTextFont(61)
        mark.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), "CMS")
        mark.SetTextSize(0.040)
        mark.SetTextFont(52)
        mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.13, 1 - (ROOT.gPad.GetTopMargin() - 0.017), " Work in Progress")
        mark.SetTextSize(0.040)
        mark.SetTextFont(42)
        mark.SetTextAlign(31)
        mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), year + " (13 TeV)")

        #ROC.SaveAs("../../fakeRateEfficiencyPlots/Run2UL_DeepCSV_Nov2022/" + year + "_ROC_WPs.pdf")
        ROC.SaveAs("Run2UL_DeepCSV_Nov2022/" + year + "_ROC_WPs.pdf")
        #ROC.Close()
       
 
if __name__ == "__main__":
    main()

