import os
import ROOT
import math
import argparse

class TaggerInfo:
 
    def __init__(self, inputFile, outputFile, wp, dataset):

        self.inputFile  = inputFile
        self.outputFile = outputFile 
        self.wp         = wp 
        self.dataset    = dataset
        self.legName = "WP_%s"%(wp)

class Eff_FakeRatePlots:

    def __init__(self):
        self.histos_    = []
        self.histoName_ = []
        self.outFile_   = []

    def defineTH1F(self, file="", branch=""):
        f = ROOT.TFile.Open(file)
        if f:
            f.cd()
            h = f.Get(branch)
            f.Close()
            return h
        else:
            print("File is nullptr")
            return 0

    def errorRatio(self, j, hnum, hdom):
        a  = hnum.GetBinContent(j)
        b  = hdom.GetBinContent(j)
        da = hnum.GetBinError(j)
        db = hdom.GetBinError(j)    

        error = 999.0
        if b != 0.0 and a != 0.0:
            error = (a/b)*math.sqrt( math.pow(da/a,2) + math.pow(db/b,2) )     

        return error
  
    def ratioTH1F(self, name, fakerateVar, Var, file, fakerate, xname, yname, join ):
        fakeRateVar = self.defineTH1F(file, fakerateVar)
        var         = self.defineTH1F(file, Var)
        fakeRateVar.Rebin(join)
        var.Rebin(join)

        fakeRate = fakeRateVar.Clone(fakerate)

        for i in range(1, fakeRate.GetNbinsX()+1):
            error = self.errorRatio(i,fakeRate,var)

            ratio = 999.0
            if var.GetBinContent(i) != 0.0:
                ratio = fakeRate.GetBinContent(i)/var.GetBinContent(i) 

            fakeRate.SetBinContent(i,ratio)
            fakeRate.SetBinError(i,error)

        fakeRate.GetXaxis().SetTitle(xname)
        fakeRate.GetYaxis().SetTitle(yname)
        fakeRate.GetYaxis().SetRangeUser(0.0, 1.2) 
        fakeRate.SetStats(0)
        self.histos_.append(fakeRate)
        self.histoName_.append("/" + fakerate)
        self.outFile_.append( "/" + name + "/") 

    def makeTH1F(self, name, rootFile, wp, dataset):
        if dataset == "TT":
            self.ratioTH1F(name, name + "/genTopMatchPtWP_" + wp + name, name + "/genTopPt_" + name, rootFile, name + "/efficiencyGenTopMatchPt", "p_{T} (GeV)", "Efficiency", 5)
        else:
            self.ratioTH1F(name, name + "/fakerateNj_"      + wp + name, name + "/nJets_"    + name, rootFile, name + "/fakerateNjets"          , "N_{jets}"   , "Fakerate"  , 1)

def makePlots(outFile, name, simpleHist, mediumHist, dataSet, simpleLeg = "WP 0.92", mediumLeg = "WP 0.95", year = "2016"):
    #Define canvas and legend
    c = ROOT.TCanvas( dataSet+name, dataSet+name, 800, 800)  
    l = ROOT.TLegend(0.30, 0.77, 0.99, 0.87) # 0.68, 0.8, 0.99, 0.9 
    topM=0.1; bottomM=0.35; rightM=0.05; leftM=0.14
    ROOT.gPad.SetAstat(0)
    ROOT.gPad.SetTopMargin(   0.0)
    ROOT.gPad.SetBottomMargin(0.0)
    ROOT.gPad.SetRightMargin( 0.0)
    ROOT.gPad.SetLeftMargin(  0.0)
    #ROOT.gPad.SetTopMargin(topM)
    #ROOT.gPad.SetBottomMargin(bottomM)
    #ROOT.gPad.SetRightMargin(rightM)
    #ROOT.gPad.SetLeftMargin(leftM)

    #Define the top and bottom TPad
    up_height     = 0.75  # please tune so that the upper figures size will meet your requirement
    dw_correction = 1.0   # 40
    dw_height     = (1.0 - up_height)*dw_correction
    pad_up = ROOT.TPad("pad_up1","pad_up1",0.0, 1.0 - up_height, 1.0,       1.0)
    pad_dw = ROOT.TPad("pad_dw1","pad_dw1",0.0, 0.0,             1.0, dw_height)
    pad_up.Draw()
    pad_dw.Draw()
  
    pad_up.SetFrameFillColor(0)
    pad_up.SetFillColor(0)
    pad_up.SetTopMargin(topM)
    pad_up.SetBottomMargin(0.02)
    pad_up.SetLeftMargin(leftM)
    pad_up.SetRightMargin(rightM)

    pad_dw.SetFrameFillColor(0)
    pad_dw.SetFillColor(0)
    pad_dw.SetTopMargin(0.04)
    pad_dw.SetBottomMargin(bottomM)
    pad_dw.SetLeftMargin(leftM)
    pad_dw.SetRightMargin(rightM)

    #Top TPad
    pad_up.cd()
    pad_up.SetGrid()
 
    absmax = max( mediumHist.GetMaximum(), simpleHist.GetMaximum() )
    mediumHist.SetMaximum( 1.2*absmax )
    mediumHist.SetMinimum(0)
    mediumHist.SetTitle("")
    mediumHist.SetTitleSize(0.002)
    mediumHist.SetTitleSize(0.05,"X")
    mediumHist.SetTitleSize(0.05,"Y")
    mediumHist.SetTitleOffset(1.2,"X")
    mediumHist.SetTitleOffset(0.5 * (0.75/0.25),"Y")
    mediumHist.SetLabelSize(0.0,"X")
    mediumHist.SetLabelSize(0.05,"Y")
    mediumHist.SetStats(False)   
    mediumHist.SetLineColor(ROOT.kCyan+2)
    mediumHist.SetLineWidth(2)
    mediumHist.Draw("hist E")

    simpleHist.SetLineColor(ROOT.kMagenta+2)
    simpleHist.SetLineWidth(2)
    simpleHist.Draw("hist E same")

    mark = ROOT.TLatex()
    mark.SetNDC(True)
    mark.SetTextAlign(11)
    mark.SetTextSize(0.06)
    mark.SetTextFont(61)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), "CMS")
    mark.SetTextSize(0.040)
    mark.SetTextFont(52)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.1, 1 - (ROOT.gPad.GetTopMargin() - 0.017), "Preliminary")
    mark.SetTextSize(0.040)
    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), (year + " (13 TeV)"))

    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.SetTextSize(0.03)    
    l.AddEntry(mediumHist, mediumLeg, "l")
    l.AddEntry(simpleHist, simpleLeg, "l")
    l.Draw()

    #Bottom TPad
    pad_dw.cd()
    pad_dw.SetGrid()

    dummy = Eff_FakeRatePlots() 
    ratio = mediumHist.Clone( (name+dataSet+"Ratio") )
    for i in range(1, ratio.GetNbinsX()+1):
        errorNum = dummy.errorRatio(i,mediumHist,simpleHist) 

        ratioNum = 999.0
        if simpleHist.GetBinContent(i) != 0.0:
            ratioNum = mediumHist.GetBinContent(i)/simpleHist.GetBinContent(i)

        ratio.SetBinContent(i,ratioNum)
        ratio.SetBinError(i,errorNum)

    maxDw = ratio.GetMaximum()
    minDw = ratio.GetMinimum()

    ratio.SetTitle("")
    ratio.SetStats(False)
    ratio.SetMaximum( 1.75 ) 
    ratio.SetMinimum( 0.25 ) 
    ratio.SetTitleOffset(1.0,"X")
    ratio.SetTitleOffset(0.3,"Y")
    ratio.SetTitleSize(0.15,"X")
    ratio.SetTitleSize(0.15,"Y")
    ratio.SetLabelSize(0.15,"X")
    ratio.SetLabelSize(0.15,"Y")
    ratio.SetNdivisions(10,"X")
    ratio.SetNdivisions( 5,"Y")
    ratio.GetYaxis().SetTitle("Ratio")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.Draw("E1")

    line = ROOT.TF1( (name+dataSet+"Line"),"1",-2000,2000)
    line.SetLineColor(ROOT.kBlue)
    line.SetLineWidth(1)
    line.Draw("same")

    outpath = ("fakeRateEfficiencyPlots/" + dataSet + name).rpartition("/")[0]
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    c.SaveAs( ( "fakeRateEfficiencyPlots/" + dataSet + name + ".pdf" ) )        

def makePlotsThesis(outFile, name, hist1, hist2, hist3, hist4, hist5, dataSet, leg1 = "WP 0.92", leg2 = "WP 0.95", leg3 = "WP 0.96", leg4 = "WP 0.97", leg5 = "WP 0.98", year = "2016"):
    #Define canvas and legend
    c = ROOT.TCanvas( (dataSet+name),(dataSet+name),800,800)  
    l = ROOT.TLegend(0.30, 0.67, 0.99, 0.87) # 0.68, 0.8, 0.99, 0.9 
    topM=0.1; bottomM=0.35; rightM=0.05; leftM=0.11
    ROOT.gPad.SetAstat(0)
    ROOT.gPad.SetTopMargin(   0.0)
    ROOT.gPad.SetBottomMargin(0.0)
    ROOT.gPad.SetRightMargin( 0.0)
    ROOT.gPad.SetLeftMargin(  0.0)
    #ROOT.gPad.SetTopMargin(topM)
    #ROOT.gPad.SetBottomMargin(bottomM)
    #ROOT.gPad.SetRightMargin(rightM)
    #ROOT.gPad.SetLeftMargin(leftM)

    #Define the top and bottom TPad
    up_height     = 1.0  # please tune so that the upper figures size will meet your requirement
    dw_correction = 1.0   # 40
    dw_height     = (1.0 - up_height)*dw_correction
    pad_up = ROOT.TPad("pad_up1","pad_up1",0.0, 1.0 - up_height, 1.0,       1.0)
    pad_dw = ROOT.TPad("pad_dw1","pad_dw1",0.0, 0.0,             1.0, dw_height)
    pad_up.Draw()
    pad_dw.Draw()
    pad_up.SetFrameFillColor(0)
    pad_up.SetFillColor(0)
    pad_up.SetTopMargin(topM)
    pad_up.SetBottomMargin(0.12)
    pad_up.SetLeftMargin(leftM)
    pad_up.SetRightMargin(rightM)
    pad_dw.SetFrameFillColor(0)
    pad_dw.SetFillColor(0)
    pad_dw.SetTopMargin(0.04)
    pad_dw.SetBottomMargin(bottomM)
    pad_dw.SetLeftMargin(leftM)
    pad_dw.SetRightMargin(rightM)

    #Top TPad
    pad_up.cd()
    pad_up.SetGrid()
    absmax = max( hist1.GetMaximum(), hist2.GetMaximum() )
    hist1.SetMaximum( 1.2*absmax )
    hist1.SetMinimum(0)
    hist1.SetTitle("")
    hist1.SetTitleSize(0.002)
    hist1.SetTitleSize(0.05,"X")
    hist1.SetTitleSize(0.05,"Y")
    hist1.SetTitleOffset(1.0,"X")
    hist1.SetTitleOffset(0.35 * (0.75/0.25),"Y")
    hist1.SetLabelSize(0.035,"X")
    hist1.SetLabelSize(0.035,"Y")
    hist1.SetStats(False)   
    hist1.SetLineColor(ROOT.TColor.GetColor("#cd5c5c"))
    hist2.SetLineColor(ROOT.TColor.GetColor("#b2e2e2"))
    hist3.SetLineColor(ROOT.TColor.GetColor("#66c2a4"))
    hist4.SetLineColor(ROOT.TColor.GetColor("#2ca25f"))
    hist5.SetLineColor(ROOT.TColor.GetColor("#00441b"))

    hist1.SetLineWidth(2)
    hist1.Draw("hist E")

    hist2.SetLineWidth(2)
    hist2.Draw("hist E same")

    hist3.SetLineWidth(2)
    hist3.Draw("hist E same")

    hist4.SetLineWidth(2)
    hist4.Draw("hist E same")

    hist5.SetLineWidth(2)
    hist5.Draw("hist E same")

    hist1.Draw("hist E same")

    mark= ROOT.TLatex()
    mark.SetNDC(True)
    mark.SetTextAlign(11)
    mark.SetTextSize(0.06)
    mark.SetTextFont(61)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), "CMS")
    mark.SetTextSize(0.040)
    mark.SetTextFont(52)
    #mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.1, 1 - (ROOT.gPad.GetTopMargin() - 0.017), "Preliminary")
    mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.1, 1 - (ROOT.gPad.GetTopMargin() - 0.017), "   Work in Progress")
    mark.SetTextSize(0.040)
    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), (year + " (13 TeV)"))

    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.SetTextSize(0.03)    
    l.AddEntry(hist1, (leg1), "l")
    l.AddEntry(hist2, (leg2), "l")
    l.AddEntry(hist3, (leg3), "l")
    l.AddEntry(hist4, (leg4), "l")
    l.AddEntry(hist5, (leg5), "l")
    l.Draw()

    outpath = ("fakeRateEfficiencyPlotsThesis/" + dataSet + name).rpartition("/")[0]
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    c.SaveAs( ( "fakeRateEfficiencyPlotsThesis/" + dataSet + name + ".pdf" ) )        

def runPlotter(taggerInfo, selections, dataSet, year):

    # -----------------------------
    # -- Looping over TaggerInfo
    # -----------------------------
    vecPlots = []
    for tI in taggerInfo:
        fakeratePlots = Eff_FakeRatePlots()
        for s in selections:
            fakeratePlots.makeTH1F(s, tI.inputFile, tI.wp, tI.dataset)

        f = ROOT.TFile(tI.outputFile,"RECREATE")
        if f.IsZombie():
            print("Cannot create " + tI.outputFile)
    
        f.cd()
    
        for i in range(0, len(fakeratePlots.histos_)):
            fakeratePlots.histos_[i].Write()
    
        f.Write()

        vecPlots.append(fakeratePlots)

    # ------------------
    # -- Making Plots
    # ------------------
    for i in range(0, len(vecPlots[0].histos_)):
        if len(vecPlots) < 3:
            makePlots(vecPlots[0].outFile_[i], vecPlots[0].histoName_[i], vecPlots[0].histos_[i], vecPlots[1].histos_[i], dataSet, taggerInfo[0].legName, taggerInfo[1].legName, year)
        else:
            makePlotsThesis(vecPlots[0].outFile_[i], vecPlots[0].histoName_[i], vecPlots[0].histos_[i], vecPlots[1].histos_[i],vecPlots[2].histos_[i], vecPlots[3].histos_[i],vecPlots[4].histos_[i], dataSet, taggerInfo[0].legName, taggerInfo[1].legName, taggerInfo[2].legName, taggerInfo[3].legName, taggerInfo[4].legName, year)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year",    dest="year",    help="Which year to run on",         type=str, required=True)
    parser.add_argument("--dataset", dest="dataset", help="Which data set to run on",     type=str, required=True)
    parser.add_argument("--inputDir",dest="inputDir",help="Input directory",              type=str, default="/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/1_ResolvedTopTagger_UL_2022/hadd_2016postVFP_ResolvedTopTagger_fakeRateEfficiency_08.04.2022")
    parser.add_argument("--wp1",     dest="wp1",     help="First working point to comp",  type=str, required=True)
    parser.add_argument("--wp2",     dest="wp2",     help="Second working point to comp", type=str, required=True)
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)
    selections = ["histos"]
    #selections = ["histos", "Njet7", "Njet8", "Njet9", "Njet10", "Njet11", "Njet12", "Njet12inc"]

    if not os.path.isdir("outputRoot"):
        os.makedirs("outputRoot")

    tagVec = []
    tagVec.append(TaggerInfo(args.inputDir + "/" + "%s_%s.root"%(args.year, args.dataset), "outputRoot/EfficiencyFakeRatePlots_%s_%s_ResolvedTopTagger.root"%(args.year, args.dataset), args.wp1, args.dataset))
    tagVec.append(TaggerInfo(args.inputDir + "/" + "%s_%s.root"%(args.year, args.dataset), "outputRoot/EfficiencyFakeRatePlots_%s_%s_ResolvedTopTagger.root"%(args.year, args.dataset), args.wp2, args.dataset))

    tagVecThesis = []
    tagVecThesis.append(TaggerInfo(args.inputDir + "/" + "%s_%s.root"%(args.year, args.dataset), "outputRoot/junk.root", "0.920", args.dataset))
    tagVecThesis.append(TaggerInfo(args.inputDir + "/" + "%s_%s.root"%(args.year, args.dataset), "outputRoot/junk.root", "0.950", args.dataset))
    tagVecThesis.append(TaggerInfo(args.inputDir + "/" + "%s_%s.root"%(args.year, args.dataset), "outputRoot/junk.root", "0.960", args.dataset))
    tagVecThesis.append(TaggerInfo(args.inputDir + "/" + "%s_%s.root"%(args.year, args.dataset), "outputRoot/junk.root", "0.970", args.dataset))
    tagVecThesis.append(TaggerInfo(args.inputDir + "/" + "%s_%s.root"%(args.year, args.dataset), "outputRoot/junk.root", "0.980", args.dataset))

    runPlotter(tagVec,       selections, "%s_%s_compareWPs_"%(args.year, args.dataset) + args.wp1 + "_" + args.wp2 + "/", args.year) 
    runPlotter(tagVecThesis, selections, "%s_%s_compareWPs_Thesis/"%(args.year, args.dataset), args.year) 

if __name__ == "__main__":
    main()
