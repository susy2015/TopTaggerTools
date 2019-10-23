#ifndef HISTOCONTAINER_H
#define HISTOCONTAINER_H

#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/MiniTupleMaker.h"

#include <vector>

#include "Math/VectorUtil.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"


class HistoContainer
{
private:
    std::vector<TH1*> histos_;
    std::string csName_;

    template<typename H, typename... Args>
    H* bookHisto(const std::string& name, Args... args)
    {
        H* hptr = new H(name.c_str(), name.c_str(), args...);
        hptr->Sumw2();
        histos_.push_back(static_cast<TH1*>(hptr));
        return hptr;
    }

public:
    TH1 *cutFlow_, *cutFlowNoWgt_, *passCuts_, *passCutsNoWgt_;

    std::vector<TH1*> jetMassNGenTopMatch;
    std::vector<TH1*> jetMassNGPMatch;

    std::vector<TH1*> jetDiscNGenTopMatch;
    std::vector<TH1*> jetDiscNGPMatch;


    HistoContainer(const std::string& csName) : csName_(csName), cutFlow_(nullptr), cutFlowNoWgt_(nullptr), passCuts_(nullptr), passCutsNoWgt_(nullptr)
    {
        for(int i = 0; i < 4; ++i)  jetMassNGenTopMatch.push_back(bookHisto<TH1D>("jetMassNGenTopMatch" + std::to_string(i),50,0,500));
        for(int i = 0; i < 6; ++i)  jetMassNGPMatch.push_back(bookHisto<TH1D>("jetMassNGPMatch" + std::to_string(i),50,0,500));

        for(int i = 0; i < 4; ++i)  jetDiscNGenTopMatch.push_back(bookHisto<TH1D>("jetDiscNGenTopMatch" + std::to_string(i),50,0,1));
        for(int i = 0; i < 6; ++i)  jetDiscNGPMatch.push_back(bookHisto<TH1D>("jetDiscNGPMatch" + std::to_string(i),50,0,1));        
    }

    ~HistoContainer()
    {
    }

    void fillWithCutFlow(const std::vector<std::pair<std::string, bool>>& cuts, const NTupleReader& tr, const float& eWeight)
    {
        const int EXTRACUT = 1; 

        if(cutFlow_ == nullptr)
        {
            cutFlow_ = bookHisto<TH1D>("cutFlow", cuts.size() + EXTRACUT, -0.5, cuts.size() + EXTRACUT - 0.5);
            cutFlowNoWgt_ = bookHisto<TH1D>("cutFlowNoWgt", cuts.size() + EXTRACUT, -0.5, cuts.size() + EXTRACUT - 0.5);
            passCuts_ = bookHisto<TH1D>("passCuts", cuts.size(), -0.5 + EXTRACUT, cuts.size() + EXTRACUT - 0.5);
            passCutsNoWgt_ = bookHisto<TH1D>("passCutsNoWgt", cuts.size(), -0.5 + EXTRACUT, cuts.size() + EXTRACUT - 0.5);

            cutFlow_->GetXaxis()->SetBinLabel(1, "all evt");
            cutFlowNoWgt_->GetXaxis()->SetBinLabel(1, "all evt");
            for(int iBin = 1; iBin <= cuts.size(); ++iBin)
            {
                cutFlow_->GetXaxis()->SetBinLabel(iBin + EXTRACUT, cuts[iBin - 1].first.c_str());
                cutFlowNoWgt_->GetXaxis()->SetBinLabel(iBin + EXTRACUT, cuts[iBin - 1].first.c_str());
                passCuts_->GetXaxis()->SetBinLabel(iBin, cuts[iBin - 1].first.c_str());
                passCutsNoWgt_->GetXaxis()->SetBinLabel(iBin, cuts[iBin - 1].first.c_str());
            }
        }

        cutFlow_->Fill(0.0, eWeight);
        cutFlowNoWgt_->Fill(0);

        bool passed = true;
        for(int i = 0; i < static_cast<int>(cuts.size()); ++i)
        {
            if(cuts[i].second)
            {
                passCuts_->Fill(i, eWeight);
                passCutsNoWgt_->Fill(i);
            }
            else
            {
                passed = false;
            }

            if(passed)
            {
                cutFlow_->Fill(i + EXTRACUT, eWeight);
                cutFlowNoWgt_->Fill(i + EXTRACUT);
            }
        }

        if(passed)
        {
            fill(tr, eWeight);
        }
    }

    void fill(const NTupleReader& tr, const float& eWeight)
    {
        const auto& fatJet_TLV   = tr.getVec_LVFromNano<float>("FatJet");
        const auto& fatJet_SDM   = tr.getVec<float>("FatJet_msoftdrop");
        const auto& fatJet_tDisc = tr.getVec<float>("FatJet_deepTag_TvsQCD");
        const auto& fatJet_nGenTopConstMatch = tr.getVec<int>("FatJet_nGenTopConstMatch");
        const auto& fatJet_nGenPartMatch = tr.getVec<int>("FatJet_nGenPartMatch");

        for(unsigned int i = 0; i < fatJet_TLV.size(); ++i)
        {
            if(fatJet_nGenTopConstMatch[i] <= 3)
            {
                jetMassNGenTopMatch[fatJet_nGenTopConstMatch[i]]->Fill(fatJet_SDM[i], eWeight);
                jetDiscNGenTopMatch[fatJet_nGenTopConstMatch[i]]->Fill(fatJet_tDisc[i], eWeight);
            }

            jetMassNGPMatch[fatJet_nGenPartMatch[i]]->Fill(fatJet_SDM[i], eWeight);
            jetDiscNGPMatch[fatJet_nGenPartMatch[i]]->Fill(fatJet_tDisc[i], eWeight);
        }
    }

    void save(TFile *f)
    {
        f->cd();
        TDirectory *td = f->mkdir(csName_.c_str(), csName_.c_str());
        td->cd();
        for(TH1* hist : histos_) hist->Write();
        f->cd();
    }
    
};

#endif
