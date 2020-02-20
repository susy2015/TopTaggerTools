from os import system, environ
import sys
from glob import glob

import datetime

sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path
from samples import SampleCollection

#sampleCollectionsToHadd = ["DYJetsToLL", "ZJetsToNuNu", "GJets", "Diboson", "TTG", "TTZ", "TTbar", "TTbarSingleLep", "WJetsToLNu", "WJetsToLNuInc", "QCD", "IncDY", "Data_JetHT", "Data_SingleMuon", "Data_SinglePhoton"]
#sampleCollectionsToHadd = ["ZJetsToNuNu", "TTbar", "TTbarSingleLep", "QCD"]
sampleCollectionsToHadd = ["TTbarSingleLep_2016", "QCD_2016"]
#sampleCollectionsToHadd = ["TTbarSingleLep",]
#sampleCollectionsToHadd = ["Data_JetHT", "Data_MET", "Data_SingleMuon", "Data_SinglePhoton"]

now = datetime.datetime.now()
directory = "submission_2020-01-22_12-22-11"

sc = SampleCollection("../sampleSets.cfg", "../sampleCollections.cfg")
scl = sc.sampleCollectionList()
for sampleCollection in scl:
    sl = sc.sampleList(sampleCollection)
    if sampleCollection in sampleCollectionsToHadd:
        files = ""
        for sample in sl:
            files += " " + " ".join(glob("%s/TT_%s_*.root"%(directory,sample[1])))
        system("hadd TT_%s-%i-%i-%i_deepResolved.root %s"%(sampleCollection, now.year, now.month, now.day, files))


