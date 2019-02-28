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


std::vector<TLorentzVector> GetHadTopLVec(const std::vector<TLorentzVector>& genDecayLVec, const std::vector<int>& genDecayPdgIdVec, const std::vector<int>& genDecayIdxVec, const std::vector<int>& genDecayMomIdxVec)
{
    std::vector<TLorentzVector> tLVec;
    for(unsigned it=0; it<genDecayLVec.size(); it++)
    {
        int pdgId = genDecayPdgIdVec.at(it);
        if(abs(pdgId)==6) //pdg(6), top quark
        {

            for(unsigned ig=0; ig<genDecayLVec.size(); ig++)
            {
                if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) ) //Is this a daughter particle of the genTop?
                {
                    int pdgId = genDecayPdgIdVec.at(ig);
                    if(abs(pdgId)==24) //pdg(24), W boson //Let's look at the W boson that is coming off of the top
                    {
                        int flag = 0;
                        for(unsigned iq=0; iq<genDecayLVec.size(); iq++)
                        {
                            if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) //Is this a daughter particle of the W?
                            {
                                int pdgid = genDecayPdgIdVec.at(iq);
                                if(abs(pdgid)== 11 || abs(pdgid)== 13 || abs(pdgid)== 15) flag++; //pdg(11), electron; pdg(13), muon; pdg(15), tau
                            }
                        }
                        if(!flag) tLVec.push_back(genDecayLVec.at(it)); //If the W didn't have a lepton daughter product then let's include in the list of hadronic tops.
                    }
                }
            }//dau. loop
        }//top cond
    }//genloop
    return tLVec;
}

std::vector<TLorentzVector> GetLepTopLVec(const std::vector<TLorentzVector>& genDecayLVec, const std::vector<int>& genDecayPdgIdVec, const std::vector<int>& genDecayIdxVec, const std::vector<int>& genDecayMomIdxVec)
{
    std::vector<TLorentzVector> tLVec;
    for(unsigned it=0; it<genDecayLVec.size(); it++)
    {
        int pdgId = genDecayPdgIdVec.at(it);
        if(abs(pdgId)==6) //pdg(6), top quark
        {
            for(unsigned ig=0; ig<genDecayLVec.size(); ig++)
            {
                if( genDecayMomIdxVec.at(ig) == genDecayIdxVec.at(it) ) //Is this a daughter particle of the genTop?
                {
                    int pdgId = genDecayPdgIdVec.at(ig);
                    if(abs(pdgId)==24) //pdg(24), W boson //Let's look at the W boson that is coming off of the top
                    {
                        int flag = 0;
                        for(unsigned iq=0; iq<genDecayLVec.size(); iq++)
                        {
                            if( genDecayMomIdxVec.at(iq) == genDecayIdxVec.at(ig) ) //Is this a daughter particle of the W?
                            {
                                int pdgid = genDecayPdgIdVec.at(iq);
                                if(abs(pdgid)== 11 || abs(pdgid)== 13 || abs(pdgid)== 15) flag++; //pdg(11), electron; pdg(13), muon; pdg(15), tau
                            }
                        }
                        if(flag) tLVec.push_back(genDecayLVec.at(it)); //If the W has a lepton daughter product then let's include it in the list of leptonic tops.
                    }
                }
            }//dau. loop
        }//top cond
    }//genloop
    return tLVec;
}

void stripRoot(std::string &path)
{
    int dot = path.rfind(".root");
    if (dot != std::string::npos)
    {
        path.resize(dot);
    }
}

float SF_13TeV(float top_pt)
{

    return exp(0.0615-0.0005*top_pt);

}

bool filterEvents(NTupleReader& tr)
{
    return true;
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

    AnaSamples::SampleSet        ss("sampleSets.cfg", runOnCondor, AnaSamples::luminosity);
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

    HistoContainer<NTupleReader> hists0Lep("Lep0"), hists1Lep("Lep1"), histsTTbar("ttbar"), histsTTbarNob("ttbarNob"), histsTTbarLep("ttbarLep"), histsQCD("QCD"), histsQCDb("QCDb"), histsLowHTQCD("lowHTQCD"), histsPhoton("photon"), histsDilepton("dilepton"), histsSimpleSemiLept("simpleSemiLep"), histsTTbar1l("histsTTbar1l"), histsTTbar2l("histsTTbar2l"), histsTTbar1lnoMET("histsTTbar1lnoMET"), histsTTbar2lnoMET("histsTTbar2lnoMET"), histsTTbarNol("histsTTbarNol");

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

            NTupleReader tr(t);

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
                const bool&  passLep20            = tr.getVar<bool>("passLep20");
                const bool&  passSingleLep30      = tr.getVar<bool>("passSingleLep30");
                const bool&  passSingleMu40       = tr.getVar<bool>("passSingleMu40");
                const bool&  passLeptonVeto       = tr.getVar<bool>("passLeptVeto");
                const bool&  passdPhis            = tr.getVar<bool>("passdPhis");
                const float& ht                   = tr.getVar<float>("HT");

                const int&    nbCSV                = tr.getVar<int>("cntCSVS");

                const bool& passMuTrigger     = tr.getVar<bool>("passMuTrigger");
                const bool& passElecTrigger   = tr.getVar<bool>("passElecTrigger");
                const bool& passMETMHTTrigger = tr.getVar<bool>("passMETMHTTrigger");
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
                    const float& puWF               = tr.getVar<float>("puWeight");
                    //std::cout << "Calculate btag WF" << std::endl;
//                    const float& bTagWF             = tr.getVar<float>((bTagSys == 1 ?  "bTagSF_EventWeightSimple_Up" : 
//                                                                       (bTagSys == -1 ? "bTagSF_EventWeightSimple_Down" : 
//                                                                                        "bTagSF_EventWeightSimple_Central")));
                    const float bTagWF = 1.0; /// FIX ME!!!!!!!!!!

                    const float& stored_weight      = tr.getVar<float>("genWeight");
                    if(enableTTbar & !noCorr)
                    {
                        const float& ttbarWF            = tr.getVar<float>("TTbarWF");
                        eWeight *= ttbarWF;
                    }
                    if(enableStored)
                    {
                        //std::cout << "weight: " << stored_weight << " jet pT: " << jetsLVec[0].Pt() << std::endl;
                        eWeight *= stored_weight;
                        //std::cout << "eWeight after stored weight: " << eWeight;
                    }
                    const float& triggerWF          = tr.getVar<float>("TriggerEffMC");

                    muTrigEff = tr.getVar<float>("muTrigWgt");

                    eWeight *= puWF * triggerWF;

                    if(!noCorr) eWeight *= bTagWF;

                    if(ttr_->getTops().size() > 0 && topReWeight){
                        for(int t = 0; t < ttr_->getTops().size(); t++){
                            eWeight *= .954;
                        }
                    }
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

                //                                    minAbsEta, maxAbsEta, minPt, maxPt
                const AnaConsts::AccRec pt45Eta24Arr = {-1,         2.4,      45,   -1  };

                int cntNJetsPt45Eta24 = AnaFunctions::countJets(jetsLVec,            pt45Eta24Arr);
                int cntNJetsPt30Eta24 = AnaFunctions::countJets(jetsLVec, AnaConsts::pt30Eta24Arr);

                const float& HT = tr.getVar<float>("HT");

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
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 1000)
                    )
                {
                    histsQCD.fill(tr, eWeight, trand);
                }

                //Low HT QCD control sample (For Fake study)
                if( (!isData || passHighHtTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 250)
                    )
                {
                    histsLowHTQCD.fill(tr, eWeight, trand);
                }

                //Simple SemiLeptonic criteria (just MC)
                if( (!isData)
                    && passNoiseEventFilter
                    && cntNJetsPt30Eta24 >= 4
                    && (ht > 250)
                    )
                {
                    histsSimpleSemiLept.fill(tr, eWeight, trand);
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
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passSingleLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && passdPhis
                    && passBLep
                    && passLepTtag
                    && deltaPhiLepMET < 0.8
                    && mTLep < 100
                    && (ht > 250)
                    && (met > 250)
                    )
                {
                    histsTTbar.fill(tr, eWeight, trand);
                }

                //TTbar1Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (met > 150)
                    )
                {
                    histsTTbar1l.fill(tr, eWeight, trand);
                }

                //TTbar2Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (met > 150)
                    )
                {
                    histsTTbar2l.fill(tr, eWeight, trand);
                }

                //TTbar No Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLeptonVeto
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && (met > 150)
                    )
                {
                    histsTTbarNol.fill(tr, eWeight, trand);
                }

                //TTbar1Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    )
                {
                    histsTTbar1lnoMET.fill(tr, eWeight, trand);
                }

                //TTbar2Lepton
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passLep20
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    )
                {
                    histsTTbar2lnoMET.fill(tr, eWeight, trand);
                }

                //semileptonic ttbar enriched control sample MET triggered
                if( (!isData || passSearchTrigger)
                    && passNoiseEventFilter
                    && passSingleLep20
                    && cntNJetsPt30Eta24 >= 4
                    && passdPhis
                    && deltaPhiLepMET < 0.8
                    && mTLep < 100
                    && (ht > 250)
                    && (met > 250)
                    )
                {
                    histsTTbarNob.fill(tr, eWeight, trand);
                }

                //semileptonic ttbar enriched control sample Mu triggered
                if( (!isData || passMuTrigger)
                    && passNoiseEventFilter
                    && passSingleMu40
                    && nbCSV >= 1
                    && cntNJetsPt30Eta24 >= 4
                    && passdPhis
                    && passBLep
                    && passLepTtag
                    && deltaPhiLepMET < 0.8
                    && mTLep < 100
                    && (ht > 200)
                    && (met > 50)
                    )
                {
                    histsTTbarLep.fill(tr, eWeight * muTrigEff, trand);
                }
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

        hists0Lep.save(f);
        hists1Lep.save(f);
        histsQCD.save(f);
        histsLowHTQCD.save(f);
        histsSimpleSemiLept.save(f);
        histsQCDb.save(f);
        histsTTbar.save(f);
        histsTTbarNob.save(f);
        histsTTbarLep.save(f);
        histsPhoton.save(f);
        histsDilepton.save(f);
        histsTTbar1l.save(f);
        histsTTbar2l.save(f);
        histsTTbarNol.save(f);
        histsTTbar1lnoMET.save(f);
        histsTTbar2lnoMET.save(f);

        f->Write();
        f->Close();
    }
}
