#!/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_8/external/slc6_amd64_gcc491/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import sys
import os
from os import system, environ
sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path

from samples import SampleCollection
import optparse 
import subprocess

filestoTransferDataMC = [environ["CMSSW_BASE"] + "/src/ResolvedTagger/Tools/DataMC", 
#                      environ["CMSSW_BASE"] + "/lib/${SCRAM_ARCH}/librecipeAUXOxbridgeMT2.so",
                      environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so", 
                      environ["CMSSW_BASE"] + "/src/ResolvedTagger/Tools/tfModel_frozen_DNN1_deepCVS_GR_balanced.pb", 
                      environ["CMSSW_BASE"] + "/src/ResolvedTagger/Tools/TopTagger.cfg",
                      environ["CMSSW_BASE"] + "/src/ResolvedTagger/Tools/sampleCollections_2017.cfg",
                      environ["CMSSW_BASE"] + "/src/ResolvedTagger/Tools/sampleSets_PostProcessed_2017.cfg"]


submitFileDataMC = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/ResolvedTagger/Tools/condorDataMC/goMakePlots.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/ResolvedTagger/Tools/condorDataMC/goMakePlots.sh,$ENV(CMSSW_BASE)/src/ResolvedTagger/Tools/condorDataMC/DataMC.tar.gz,$ENV(CMSSW_BASE)/src/ResolvedTagger/Tools/condorDataMC/$ENV(CMSSW_VERSION).tar.gz
Output = logs/makePlots_$(Process).stdout
Error = logs/makePlots_$(Process).stderr
Log = logs/makePlots_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)
"""

parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-n', dest='numfile', type='int', default = 5, help="number of files per job")
parser.add_option ('-d', dest='datasets', type='string', default = '', help="List of datasets 'Signal_T1tttt_mGluino2000_mLSP100'")
parser.add_option ('-l', dest='dataCollections', action='store_true', default = False, help="List all datacollections")
parser.add_option ('-L', dest='dataCollectionslong', action='store_true', default = False, help="List all datacollections and sub collections")
#parser.add_option ('-f', dest='subsample', type='string', default = '', help="Nmae of subsample 'WJetsToLNu_HT_2500toInf'")

options, args = parser.parse_args()

submitFile = ""
exeName = ""

def makeExeAndFriendsTarrball(filestoTransfer, fname):
    if not options.dataCollections and not options.dataCollectionslong:
        #WORLDSWORSESOLUTIONTOAPROBLEM
        system("mkdir -p WORLDSWORSESOLUTIONTOAPROBLEM")
        for fn in filestoTransfer:
            system("cd WORLDSWORSESOLUTIONTOAPROBLEM; ln -s %s"%fn)
        
        tarallinputs = "tar czvf %s.tar.gz WORLDSWORSESOLUTIONTOAPROBLEM --dereference"%fname
        print tarallinputs
        system(tarallinputs)
        system("rm -r WORLDSWORSESOLUTIONTOAPROBLEM")

if not options.dataCollections and not options.dataCollectionslong:
    system("tar --exclude-caches-all --exclude-vcs -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp")

exeName = "DataMC"
submitFile = submitFileDataMC
makeExeAndFriendsTarrball(filestoTransferDataMC, "DataMC")


nFilesPerJob = options.numfile
#subsamplename = options.subsample

fileParts = [submitFile]
sc = SampleCollection("../sampleSets_PostProcessed_2017.cfg", "../sampleCollections_2017.cfg")
datasets = []

if options.datasets:
    datasets = options.datasets.split(',')
else:
    print "No dataset pecified"
    exit(0)

for ds in datasets:
    ds = ds.strip()

    for s, n, e in sc.sampleList(ds):
        print n
        print s
       # if n!= subsamplename:
           # continue
        f = open(s)
        if not f == None:
            count = 0
            for l in f:
                if '.root' in l and not 'failed' in l:
                    count = count + 1
            for startFileNum in xrange(0, count, nFilesPerJob):
                #fileParts.append("Arguments = %s %s %i %i %s\nQueue\n\n"%(n, rel_name, startFileNum, nFilesPerJob, s))
                fileParts.append("Arguments = %s $ENV(CMSSW_VERSION) %i %i %s\nQueue\n\n"%(n, startFileNum, nFilesPerJob, s))
            f.close()

fout = open("condorDataMC_submit.txt", "w")
fout.write(''.join(fileParts))
fout.close()

system('mkdir -p logs')
system("echo 'condor_submit condorDataMC_submit.txt'")
system('condor_submit condorDataMC_submit.txt')
