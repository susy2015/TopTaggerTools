#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/SATException.h"

#include "TopTaggerTools/Tools/include/HistoContainer.h"

#include "derivedTupleVariables.h"
#include "SusyAnaTools/Tools/BTagCorrector.h"
#include "SusyAnaTools/Tools/TTbarCorrector.h"
#include "SusyAnaTools/Tools/ISRCorrector.h"
#include "SusyAnaTools/Tools/PileupWeights.h"
#include "SusyAnaTools/Tools/customize.h"

#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/Constituent.h"

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
    AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);

    if(dataSets.find("Data") != std::string::npos)
    {
       doWgt = false;
    }

    if(dataSets.find("TT") != std::string::npos)
    {
       enableTTbar = true;
    }

    if(dataSets.find("Pt15to7000") != std::string::npos && !overrideStored)
    {
       enableStored = true;
    }

    HistoContainer<NTupleReader> histsTTbar("ttbar"), histsTTbarLep("ttbarLep"), histsQCD("QCD"), histsQCDb("QCDb"), histsLowHTQCD("lowHTQCD"), histsLowHTQCDb("lowHTQCDb"), histsPhoton("photon"), histsDilepton("dilepton");

    TRandom* trand = new TRandom3();

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
            plotterFunctions::PrepareTopCRSelection prepTopCR(JECSys);
            plotterFunctions::PrepareTopVars prepareTopVars("TopTagger.cfg");
            plotterFunctions::TriggerInfo triggerInfo(false, false);
            TTbarCorrector ttbarCorrector(false, "");
            ISRCorrector ISRcorrector("allINone_ISRJets.root","","");

            NTupleReader tr(t, {"run"});

            tr.registerFunction(filterEvents);
            tr.registerFunction(prepTopCR);
            tr.registerFunction(prepareTopVars);
            tr.registerFunction(triggerInfo);
            tr.registerFunction(ttbarCorrector);
            tr.registerFunction(ISRcorrector);

            float fileWgt = fs.getWeight();

            const int printInterval = 1000;
            int printNumber = 0;

            while(tr.getNextEvent())
            {
                if(nEvts > 0 && tr.getEvtNum() > nEvts) break;
                if(tr.getEvtNum() / printInterval > printNumber)
                {
                    printNumber = tr.getEvtNum() / printInterval;
                    std::cout << "Event #: " << printNumber * printInterval << std::endl;
                }

                const float& met    = tr.getVar<float>("MET_pt");
                const float& metphi = tr.getVar<float>("MET_phi");

                TLorentzVector MET;
                MET.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

                const bool&  passNoiseEventFilter = (tr.getVar<bool>("passNoiseEventFilter") || noNoiseFilter);
                const bool&  passSingleLep20      = tr.getVar<bool>("passSingleLep20");
                const bool&  passSingleMu40       = tr.getVar<bool>("passSingleMu40");
                const bool&  passLeptonVeto       = tr.getVar<bool>("passLeptVeto");
                const bool&  passdPhis            = tr.getVar<bool>("passdPhis");
                const float& ht                   = tr.getVar<float>("HT");

                const int&    nbCSV                = tr.getVar<int>("cntCSVS");

                const bool& passMuTrigger     = tr.getVar<bool>("passMuTrigger");
                const bool& passElecTrigger   = tr.getVar<bool>("passElecTrigger");
                const bool& passSearchTrigger = tr.getVar<bool>("passSearchTrigger");
                const bool& passHighHtTrigger = tr.getVar<bool>("passHighHtTrigger");
                const bool& passPhotonTrigger = tr.getVar<bool>("passPhotonTrigger");

                const bool& passfloatLep      = tr.getVar<bool>("passfloatLep");

                const std::vector<TLorentzVector>& cutMuVec = tr.getVec<TLorentzVector>("cutMuVec");
                const std::vector<float>& cutMuMTlepVec = tr.getVec<float>("cutMuMTlepVec");
                const std::vector<TLorentzVector>& cutElecVec = tr.getVec<TLorentzVector>("cutElecVec");
                const std::vector<float>& cutElecMTlepVec = tr.getVec<float>("cutElecMTlepVec");

                const TopTaggerResults *ttr_                 =  tr.getVar<TopTaggerResults const*>("ttrMVA");

                const float isData = !tr.checkBranch("GenPart_pt");

                float eWeight = fileWgt;

                float muTrigEff = 1.0;
                const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>(jetVecLabel);

                if(!isData && doWgt)
                {
                    const float& puWF               = 1.0;//tr.getVar<float>("puWeight");
                    //std::cout << "Calculate btag WF" << std::endl;
//                    const float& bTagWF             = tr.getVar<float>((bTagSys == 1 ?  "bTagSF_EventWeightSimple_Up" : 
//                                                                       (bTagSys == -1 ? "bTagSF_EventWeightSimple_Down" : 
//                                                                                        "bTagSF_EventWeightSimple_Central")));
                    const float bTagWF = 1.0; /// FIX ME!!!!!!!!!!

                    const float& stored_weight      = tr.getVar<float>("genWeight_sign");
                    if(enableTTbar & !noCorr)
                    {
                        const float& ttbarWF            = 1.0;//tr.getVar<float>("TTbarWF");
                        eWeight *= ttbarWF;
                    }
//                    const float& triggerWF          = tr.getVar<float>("TriggerEffMC");

                    muTrigEff = tr.getVar<float>("muTrigWgt");

                    eWeight *= stored_weight;
                    eWeight *= puWF;
                    eWeight *= bTagWF;
                }
                
                const std::vector<float>& recoJetsBtag     = tr.getVec<float>("Jet_btagDeepB");

                //Find lepton (here it is assumed there is exactly 1 lepton)
                TLorentzVector lepton;
                float mTLep = 999.9;
                for(unsigned int i = 0; i < cutMuVec.size(); ++i)
                {
                    if(cutMuVec[i].Pt() > 20)
                    {
                        lepton = cutMuVec[i];
                        mTLep = cutMuMTlepVec[i];
                        break;
                    }
                }
                for(unsigned int i = 0; i < cutElecVec.size(); ++i)
                {
                    if(cutElecVec[i].Pt() > 20)
                    {
                        lepton = cutElecVec[i];
                        mTLep = cutElecMTlepVec[i];
                        break;
                    }
                }

                tr.registerDerivedVar("lepton", lepton);
                tr.registerDerivedVar("mTLep", mTLep);

                int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);

                const bool& passPhoton200 = tr.getVar<bool>("passPhoton200");

                // calculate passBLep
                bool passBLep = false;
                bool passLepTtag = false;
                for(int i = 0; i < jetsLVec.size(); i++)
                {
                    //Is this a b-tagged jet (loose wp?)?
                    if(recoJetsBtag[i] < 0.8) continue;

                    float lepTopMass = (lepton + jetsLVec[i]).M();
                    if(lepTopMass > 30 && lepTopMass < 180)
                    {
                        passLepTtag = true;
                    }

                    if(jetsLVec[i].DeltaR(lepton) < 1.5)
                    {
                        passBLep = true;
                    }
                }

                float deltaPhiLepMET = fabs(lepton.DeltaPhi(MET));

                //High HT QCD control sample
                std::vector<std::pair<std::string, bool>> htQCDCuts = {
                    {"trig",    (!isData || passHighHtTrigger)},
                    {"filter",  passNoiseEventFilter},
                    {"lepVeto", passLeptonVeto},
                    {"nJet",    cntNJetsPt30Eta24 >= 4},
                    {"HT1000",  ht > 1000},
                };
                histsQCD.fillWithCutFlow(htQCDCuts, tr, eWeight, trand);

                //Low HT QCD control sample (For Fake study)
                std::vector<std::pair<std::string, bool>> lowHTQCDCuts = {
                    {"trig",    (!isData || passHighHtTrigger)},
                    {"filter",  passNoiseEventFilter},
                    {"lepVeto", passLeptonVeto},
                    {"nJet",    cntNJetsPt30Eta24 >= 4},
                    {"HT250",   ht > 250},
                };
                histsLowHTQCD.fillWithCutFlow(lowHTQCDCuts, tr, eWeight, trand);


                //Low HT QCD control sample (For Fake study)
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 250)
                    )
                {
                    histsLowHTQCDb.fill(tr, eWeight, trand);
                }

                //High HT QCD + b control sample
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 1000)
                    )
                {
                    histsQCDb.fill(tr, eWeight, trand);
                }

                //photon control sample
                if( (!isData || passPhotonTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && cntNJetsPt30Eta24 >= 4
                    && passPhoton200
                    && ht > 400
                    )
                {
                    histsPhoton.fill(tr, eWeight, trand);
                }

                //dilepton control sample
                if( (!isData || passMuTrigger)
                    && passNoiseEventFilter
                    && cntNJetsPt30Eta24 >= 4
                    && passfloatLep                    
                    )
                {
                    histsDilepton.fill(tr, eWeight * muTrigEff, trand);
                }

                //semileptonic ttbar enriched control sample MET triggered
                std::vector<std::pair<std::string, bool>> ttbarCuts = {
                    {"trig",     (!isData || passSearchTrigger)},
                    {"filter",   passNoiseEventFilter},
                    {"njets",    passSingleLep20},
                    {"mu40",     nbCSV >= 1},
                    {"mTLep",    cntNJetsPt30Eta24 >= 4},
                    {"nb",       passdPhis},
                    {"dphi",     passBLep},
                    {"BLep",     passLepTtag},
                    {"LepTTag",  deltaPhiLepMET < 0.8},
                    {"dPhiLMET", mTLep < 100},
                    {"HT250",    ht > 250},
                    {"MET250",   met > 250}
                };
                histsTTbar.fillWithCutFlow(ttbarCuts, tr, eWeight, trand);

                //semileptonic ttbar enriched control sample Mu triggered
                std::vector<std::pair<std::string, bool>> ttbarLepCuts = {
                    {"trig",     (!isData || passMuTrigger)},
                    {"filter",   passNoiseEventFilter},
                    {"njets",    cntNJetsPt30Eta24 >= 4},
                    {"mu40",     passSingleMu40},
                    {"mTLep",    mTLep < 100},
                    {"nb",       nbCSV >= 1},
                    {"dphi",     passdPhis},
                    {"BLep",     passBLep},
                    {"LepTTag",  passLepTtag},
                    {"dPhiLMET", deltaPhiLepMET < 0.8},
                    {"HT200",    ht > 200},
                    {"MET50",    met > 50} 
                };
                histsTTbarLep.fillWithCutFlow(ttbarLepCuts, tr, eWeight * muTrigEff, trand);
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

        histsQCD.save(f);
        histsLowHTQCD.save(f);
        histsLowHTQCDb.save(f);
        histsQCDb.save(f);
        histsTTbar.save(f);
        histsTTbarLep.save(f);
        histsPhoton.save(f);
        histsDilepton.save(f);

        f->Write();
        f->Close();
    }
}
