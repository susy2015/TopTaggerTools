#include "plot.h"


int main()
{   
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry  root file                 draw options  draw color
    histInfo data = {"Data",    "Data_SingleMuon-combined.root", "PEX0",       kBlack};

    //vector summarizing background histograms to include in the plot
    //vector summarizing background histograms to include in the plot 
    std::vector<histInfo> bgEntries = {
//        {"GJets",         "GJets-combined.root",      "hist", kViolet},
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},   
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
//        {"DYJets",        "DYJets-combined.root",     "hist", kCyan},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"tW",            "tW-combined.root",         "hist", kSpring}, 
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta}
    };

    //These are the background that we want to subtract from the data
    std::vector<histInfo> unwantedEntries = {
        {"QCD",           "QCD-combined.root",        "hist", kBlue+2},
        {"Rare",          "Rare-combined.root",       "hist", kAzure},
        {"WJets",         "WJets-combined.root",      "hist", kTeal},
        {"Diboson",       "Diboson-combined.root",    "hist", kPink},
        {"tW",            "tW-combined.root",         "hist", kSpring}, 
        {"ZJets",         "ZJets-combined.root",      "hist", kMagenta+3}
    };

    //These are the entires that have the objects we want (hadronically decaying tops)
    std::vector<histInfo> wantedEntries = {
        {"TTbar",         "TTbar-combined.root",      "hist", kOrange},
        {"TTZ",           "TTZ-combined.root",        "hist", kMagenta}
    };

    //vector summarizing signal histograms to include in the plot
//    std::vector<histInfo> sigEntries = {
//        {"T2tt (1000, 1)", "myhistos/Signal_fastsim_T2tt_mStop-1000.root", "hist", kGreen + 2},
//    };

    //make plotter object with the required sources for histograms specified
//    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));
    Plotter plt(std::move(data), std::move(bgEntries));
    Plotter pltSub(std::move(data), std::move(unwantedEntries), std::move(wantedEntries)); //Plotter to make background subtracted plots

    plt.plot("ttbarLep/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/MET", "p_{T}^{miss} [GeV]", "Events", true, -1, -1, 5);

    plt.plot("ttbarLep/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/nJets", "Jets/Event", "Events", true);
    plt.plot("ttbarLep/nVertices", "primary vertices/event", "Events", true);
    plt.plot("ttbarLep/photon", "Photon p_{T} [GeV]", "Events", true);
    plt.plot("ttbarLep/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);

    pltSub.plot("ttbarLep/HT", "H_{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/nJets", "Jets/Event", "Events", true);
    pltSub.plot("ttbarLep/nVertices", "primary vertices/event", "Events", true);
    pltSub.plot("ttbarLep/photon", "Photon p_{T} [GeV]", "Events", true);
    pltSub.plot("ttbarLep/topMass", "m_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/topP", "p_{top} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/topPt", "p_{top}^{T} [GeV]", "Events", true, -1, -1, 5);
//    plt.plot("ttbarLep/DiTopMass", "Invariant mass of ditop events [GeV]", "Events", true, -1, -1, 5);

    plt.plot("ttbarLep/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    plt.plot("ttbarLep/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);

    pltSub.plot("ttbarLep/bestTopPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/bestTopMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/bestTopEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/bestTopCandMass", "m_{top,best} [GeV] (passes top tagger)", "Events", true, -1, -1, 5);
    pltSub.plot("ttbarLep/bestTopCandEta", "#eta_{top,best} (passes top tagger)", "Events", true, -1, -1, 5);

    plt.plotEff(true, "ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    pltSub.plotEff(true, "ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Events", true, -1, -1, 5);

    plt.plotEff(false, "ttbarLep/genTopMatchPt", "ttbarLep/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/genTopMatchMass", "ttbarLep/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/genTopMatchEta", "ttbarLep/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

    pltSub.plotEff(false, "ttbarLep/genTopMatchPt", "ttbarLep/genTopPt", "p_{genTop}^{T} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/genTopMatchMass", "ttbarLep/genTopMass", "m_{genTop} [GeV]", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/genTopMatchEta", "ttbarLep/genTopEta", "#eta_{genTop}", "Events", true, -1, -1, 5);

    plt.plotSF("ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    pltSub.plotSF("ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Tagging Rate", false, -1, -1, 5);

    plt.plotOnlySF("ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    pltSub.plotOnlySF("ttbarLep/bestTopPt", "ttbarLep/bestTopCandPt", "p_{top,best}^{T} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbarLep/bestTopMass", "ttbarLep/bestTopCandMass", "m_{top,best} [GeV]", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbarLep/bestTopEta", "ttbarLep/bestTopCandEta", "#eta_{top,best}", "Scale Factor", false, -1, -1, 5);

    plt.plotEff(true, "ttbarLep/METTagged",        "ttbarLep/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbarLep/HTTagged",         "ttbarLep/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbarLep/nJetsTagged",      "ttbarLep/nJets",     "# of Jets per Event",     "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbarLep/nVerticesTagged",  "ttbarLep/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(true, "ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotEff(false, "ttbarLep/METTagged",        "ttbarLep/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/HTTagged",         "ttbarLep/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/nJetsTagged",      "ttbarLep/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    plt.plotEff(false, "ttbarLep/nVerticesTagged",  "ttbarLep/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    plt.plotEff(false, "ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(true, "ttbarLep/METTagged",        "ttbarLep/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbarLep/HTTagged",         "ttbarLep/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbarLep/nJetsTagged",      "ttbarLep/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(true, "ttbarLep/nVerticesTagged",  "ttbarLep/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(true, "ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    pltSub.plotEff(false, "ttbarLep/METTagged",        "ttbarLep/MET",       "p_{miss}^{T} [GeV]",      "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/HTTagged",         "ttbarLep/HT",        "H^{T} [GeV]",             "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/nJetsTagged",      "ttbarLep/nJets",     "# of Jets per Event",     "Events", true, -1, -1, -1);
    pltSub.plotEff(false, "ttbarLep/nVerticesTagged",  "ttbarLep/nVertices", "# of Vertices per Event", "Events", true, -1, -1, 5);
    pltSub.plotEff(false, "ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Events", true, -1, -1, 5);
    
    plt.plotSF("ttbarLep/METTagged",       "ttbarLep/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbarLep/HTTagged",        "ttbarLep/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbarLep/nJetsTagged",     "ttbarLep/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    plt.plotSF("ttbarLep/nVerticesTagged", "ttbarLep/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    plt.plotSF("ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    pltSub.plotSF("ttbarLep/METTagged",       "ttbarLep/MET",       "p_{miss}^{T} (GeV)",      "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbarLep/HTTagged",        "ttbarLep/HT",        "H^{T} (GeV)",             "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbarLep/nJetsTagged",     "ttbarLep/nJets",     "# of Jets per Event",     "Tagging Rate", false, -1, -1, -1);
    pltSub.plotSF("ttbarLep/nVerticesTagged", "ttbarLep/nVertices", "# of Vertices per Event", "Tagging Rate", false, -1, -1, 5);
    pltSub.plotSF("ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Tagging Rate", false, -1, 1, 5);

    plt.plotOnlySF("ttbarLep/METTagged",       "ttbarLep/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbarLep/HTTagged",        "ttbarLep/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbarLep/nJetsTagged",     "ttbarLep/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    plt.plotOnlySF("ttbarLep/nVerticesTagged", "ttbarLep/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    plt.plotOnlySF("ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

    pltSub.plotOnlySF("ttbarLep/METTagged",       "ttbarLep/MET",       "p_{miss}^{T} (GeV)",      "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbarLep/HTTagged",        "ttbarLep/HT",        "H^{T} (GeV)",             "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbarLep/nJetsTagged",     "ttbarLep/nJets",     "# of Jets per Event",     "Scale Factor", false, -1, -1, -1);
    pltSub.plotOnlySF("ttbarLep/nVerticesTagged", "ttbarLep/nVertices", "# of Vertices per Event", "Scale Factor", false, -1, -1, 5);
    pltSub.plotOnlySF("ttbarLep/photonTagged",     "ttbarLep/photon",    "Photon p_{T} [GeV]", "Scale Factor", false, -1, 1, 5);

}
