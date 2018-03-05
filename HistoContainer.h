#include "../../SusyAnaTools/Tools/NTupleReader.h"
#include "../../SusyAnaTools/Tools/samples.h"
#include "../../SusyAnaTools/Tools/SATException.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "BTagCorrector.h"
#include "TTbarCorrector.h"
#include "ISRCorrector.h"
#include "PileupWeights.h"
#include "customize.h"

#include "TopTaggerResults.h"
#include "Constituent.h"

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include "math.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"

class HistoContainer
{
private:
    std::vector<TH1*> histos;
    std::string csName_;

    template<typename H, typename... Args>
    H* bookHisto(const std::string& name, Args... args)
    {
        H* hptr = new H((csName_ + name).c_str(), (csName_ + name).c_str(), args...);
        hptr->Sumw2();
        histos.push_back(static_cast<TH1*>(hptr));
        return hptr;
    }

public:
    TH1* hMET;
    TH1* hNJets;
    TH1* hNBJets;
    TH1* hNVertices;
    TH1* hTopMass;
    TH1* hTopP;
    TH1* hTopPt;
    TH1* hDiTopMass;
    TH1 *topPt, *topMass, *topEta;
    TH1 *topCandPt, *topCandMass, *topCandEta;
    TH1 *genTopPt, *genTopMass, *genTopEta;
    TH1 *genTopMatchPt, *genTopMatchMass, *genTopMatchEta;
    TH1 *bestTopCandPt, *bestTopCandMass, *bestTopCandEta;
    TH1 *bestTopCandSumPt, *bestTopCandSumMass, *bestTopCandSumEta;
    TH1 *bestTopGenPt, *bestTopGenMass, *bestTopGenEta;
    TH1 *bestTopNotGenPt, *bestTopNotGenMass, *bestTopNotGenEta;
    TH1 *bestTopPt, *bestTopMass, *bestTopEta;
    TH1 *randomTopCandPt,   *randomTopCandMass,   *randomTopCandEta;
    TH1 *randomTopPt, *randomTopMass, *randomTopEta;
    TH1 *fakerateMET, *fakerateNj, *fakerateNb;
    TH1 *fakerateMET2, *fakerateNj2, *fakerateNb2, *fakerateNvert2;

    TH1 *massTemplateTop, *massTemplateNotTop;

    TH1 *allSumPt, *bestSumPt, *genSumPt;

    TH2 *topCandMassByPt, *massTemplateTopByPt, *massTemplateNotTopByPt;
    TH2 *bestTopCandSumMassByPt;
    TH2 *bestTopCandSumMassRecoMatchByPt;
    TH2 *massTemplateGen0MatchByPt;
    TH2 *massTemplateGen1MatchByPt;
    TH2 *massTemplateGen2MatchByPt;
    TH2 *massTemplateGen3MatchByPt;

    HistoContainer(const std::string& csName);
    void fill(const NTupleReader& tr, const double& eWeight, TRandom* trand);
    void save(TFile *f);

