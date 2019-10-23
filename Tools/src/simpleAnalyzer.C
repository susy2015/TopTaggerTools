#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/SATException.h"

#include "TopTaggerTools/Tools/include/HistoContainer.h"

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include "math.h"

#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TChain.h"

void stripRoot(std::string &path)
{
    int dot = path.rfind(".root");
    if (dot != std::string::npos)
    {
        path.resize(dot);
    }
}

bool filterEvents(NTupleReader& tr)
{
    auto& jet_pt = tr.getVec<float>("Jet_pt");
    return jet_pt.size() >= 4 && jet_pt[3] >= 20;
}

int main(int argc, char* argv[])
{

    std::string jetVecLabel           = "jetsLVec";

    int opt;
    int option_index = 0;

    static struct option long_options[] = {
        {"condor",             no_argument, 0, 'c'},
        {"TTbar weight",       no_argument, 0, 't'},
        {"Stored weight",      no_argument, 0, 's'},
        {"Supress Noise Filter", no_argument, 0, 'f'},
        {"Override Stored weight", no_argument, 0, 'r'},
        {"No TTbar or Btag corrections", no_argument, 0, 'a'},
        {"Apply tight TopTag reweighting", no_argument, 0, 'R'},
        {"Use prodjetNolep branches",no_argument, 0, 'z'},
        {"no event weighting", no_argument, 0, 'd'},
        {"jecUnc up",          no_argument, 0, 'u'},
        {"jecUnc down",        no_argument, 0, 'b'}, 
        {"bTagUnc up",         no_argument, 0, 'U'},
        {"bTagUnc down",       no_argument, 0, 'B'}, 
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
        {"numEvts",      required_argument, 0, 'T'},
        {"output",       required_argument, 0, 'O'}
    };

    bool runOnCondor = false, enableTTbar = false, doWgt = true, enableStored = false, overrideStored = false, noCorr = false, useAltBranch = false, noNoiseFilter = false, onlyTTbarAllHad = false, onlyTTbarSingleLep = false, onlyTTbarDiLep = false, topReWeight = false;
    int nFiles = -1, startFile = 0, nEvts = -1, tF = -1;
    std::string dataSets = "Signal_T2tt_mStop850_mLSP100", filename = "example.root";

    int JECSys = 0, bTagSys = 0;

    while((opt = getopt_long(argc, argv, "ctsfraRzdubUBD:N:M:E:O:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'c':
            runOnCondor = true;
            break;

        case 't':
            enableTTbar = true;
            break;

        case 's':
            enableStored = true;
            break;

        case 'f':
            noNoiseFilter = true;
            break;

        case 'r':
            overrideStored = true;
            enableStored = false;
            break;

        case 'a':
            noCorr = true;
            break;

        case 'R':
            topReWeight = true;
            break;

       case 'z':
            useAltBranch = true;
            break;

        case 'd':
            doWgt = false;
            break;

        case 'u':
            JECSys = 1;
            break;

        case 'b':
            JECSys = -1;
            break;

        case 'U':
            bTagSys = 1;
            break;

        case 'B':
            bTagSys = -1;
            break;

        case 'D':
            dataSets = optarg;
            break;

        case 'N':
            nFiles = int(atoi(optarg));
            break;

        case 'M':
            startFile = int(atoi(optarg));
            break;

        case 'E':
            nEvts = int(atoi(optarg));
            break;

        case 'O':
            filename = optarg;
            std::cout << "Filename: " << filename << std::endl;

        }
    }

    if(JECSys == 1) std::cout << "JEC uncertainty up." << std::endl;
    if(JECSys == -1) std::cout << "JEC uncertainty down." << std::endl;

    if(bTagSys == 1) std::cout << "bTag uncertainty up." << std::endl;
    if(bTagSys == -1) std::cout << "bTag uncertainty down." << std::endl;

    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        stripRoot(filename);
        sprintf(thistFile, "%s_%s_%d.root", filename.c_str(), dataSets.c_str(), startFile);
        filename = thistFile;
    }

    TH1::AddDirectory(false);

    bool savefile = true;
    if(filename == "-")
    {
        savefile = false;
    }

    AnaSamples::SampleSet        ss("sampleSets.cfg", runOnCondor);
    //AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);

    HistoContainer histsBaselineHighDM("baselineHighDm");
    HistoContainer histsInclusive("Inclusive");

    MiniTupleMaker mtm("miniTree_" + filename, "Events");
    mtm.setTupleVars({"Stop0l_evtWeight", "Pass_JetID", "Pass_CaloMETRatio", "Pass_EventFilter", "Pass_highDM", "Pass_lowDM", "Pass_Baseline", "FatJet_pt", "FatJet_eta", "FatJet_phi", "FatJet_mass", "FatJet_msoftdrop", "FatJet_deepTag_TvsQCD", "FatJet_nGenPartMatch", "FatJet_nGenTopConstMatch", "FatJet_Stop0l"});

    try
    {

        //for(auto& fs : sc[dataSets])
        auto& fs = ss[dataSets];
        {
            TChain *t = new TChain(fs.treePath.c_str());
            fs.addFilesToChain(t, startFile, nFiles);

            std::cout << "File: " << fs.filePath << std::endl;
            std::cout << "Tree: " << fs.treePath << std::endl;

            //plotterFunctions::SystematicPrep sysPrep;
            NTupleReader tr(t, {"run"});

            float fileWgt = fs.getWeight();

            const int printInterval = 10000;
            int printNumber = 0;


            while(tr.getNextEvent())
            {
                if(nEvts > 0 && tr.getEvtNum() > nEvts) break;
                if(tr.getEvtNum() / printInterval > printNumber)
                {
                    printNumber = tr.getEvtNum() / printInterval;
                    std::cout << "Event #: " << printNumber * printInterval << std::endl;
                }

                double eWeight = fileWgt;
                
                //cuts
                const bool& Pass_JetID        = tr.getVar<bool>("Pass_JetID");
                const bool& Pass_CaloMETRatio = tr.getVar<bool>("Pass_CaloMETRatio");
                const bool& Pass_EventFilter  = tr.getVar<bool>("Pass_EventFilter");
                const bool& Pass_highDM       = tr.getVar<bool>("Pass_highDM");

                const auto& fatJet_TLV   = tr.getVec_LVFromNano<float>("FatJet");

                const auto& genPart_TLV              = tr.getVec_LVFromNano<float>("GenPart");
                const auto& genPart_pdgId            = tr.getVec<int>("GenPart_pdgId");
                const auto& genPart_genPartIdxMother = tr.getVec<int>("GenPart_genPartIdxMother");
                const auto& genPart_statusFlags      = tr.getVec<int>("GenPart_statusFlags");

                std::pair<std::vector<TLorentzVector>, std::vector<std::vector<const TLorentzVector*>>> genTops = ttUtility::GetTopdauGenLVecFromNano(genPart_TLV, genPart_pdgId, genPart_statusFlags, genPart_genPartIdxMother);

                std::vector<int>& fatJet_nGenTopConstMatch = tr.createDerivedVec<int>("FatJet_nGenTopConstMatch", fatJet_TLV.size());
                for(unsigned int i = 0; i < fatJet_TLV.size(); ++i)
                {
                    for(const auto& genTopConst : genTops.second)
                    {
                        int nMatch = 0;
                        for(const auto& gp : genTopConst)
                        {
                            if(ROOT::Math::VectorUtil::DeltaR(fatJet_TLV[i], *gp) < 0.6)
                            {
                                ++nMatch;
                            }
                        }
                        if(nMatch > fatJet_nGenTopConstMatch[i]) fatJet_nGenTopConstMatch[i] = nMatch;
                    }
                }

                std::vector<int>& fatJet_nGenPartMatch = tr.createDerivedVec<int>("FatJet_nGenPartMatch", fatJet_TLV.size());
                for(unsigned int i = 0; i < fatJet_TLV.size(); ++i)
                {
                    for(unsigned int j = 0; j < genPart_TLV.size(); ++j)
                    {
                        const int absPdgId = abs(genPart_pdgId[j]);
                        if((absPdgId >= 1 && absPdgId <= 5) || absPdgId == 21) 
                        {
                            if(ROOT::Math::VectorUtil::DeltaR(fatJet_TLV[i], genPart_TLV[j]) < 0.6)
                            {
                                ++fatJet_nGenPartMatch[i];
                            }
                        }
                    }
                    if(fatJet_nGenPartMatch[i] >= 6) fatJet_nGenPartMatch[i] = 5;
                }

                //High DM region
                std::vector<std::pair<std::string, bool>> highDMCuts = {
                    {"JetID",        Pass_JetID},
                    {"CaloMETRatio", Pass_CaloMETRatio},
                    {"EventFilter",  Pass_EventFilter},
                    {"highDM",       Pass_highDM},
                };
                histsBaselineHighDM.fillWithCutFlow(highDMCuts, tr, eWeight);

                //Inclusive
                std::vector<std::pair<std::string, bool>> eventFilterOnly = {
                    {"JetID",        Pass_JetID},
                    {"CaloMETRatio", Pass_CaloMETRatio},
                    {"EventFilter",  Pass_EventFilter},
                };
                histsInclusive.fillWithCutFlow(eventFilterOnly, tr, eWeight);

                if(tr.isFirstEvent())
                {
                    mtm.initBranches(tr);
                }

                if(tr.getVar<bool>("Pass_Baseline")) mtm.fill();
            }
        }
    }
    catch(const std::string e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const TTException e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const SATException e)
    {
        std::cout << e << std::endl;
        return 0;
    }


    if(savefile)
    {
        std::cout << "Saving root file..." << std::endl;

        TFile *f = new TFile(filename.c_str(),"RECREATE");
        if(f->IsZombie())
        {
            std::cout << "Cannot create " << filename << std::endl;
            throw "File is zombie";
        }

        histsBaselineHighDM.save(f);
        histsInclusive.save(f);

        f->Write();
        f->Close();
    }
}
