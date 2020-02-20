#!/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_8/external/slc6_amd64_gcc491/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import sys
import os
from os import system, environ
sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path

from samples import SampleCollection
import optparse 
import subprocess
import datetime

# TopTagger.cfg
mvaFileName = ""
with file(environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/TopTagger.cfg") as meowttcfgFile:
    for line in meowttcfgFile:
        line = line.split("#")[0]
        if "modelFile" in line:
            mvaFileName = line.split("=")[1].strip().strip("\"")
            break

#here I hack in the tarball for GMP, this needs to be generalized to the other options 

#Here is the configuration for the Data/MC validation of the TopTagger 
filestoTransferTT  = [environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/simpleAnalyzer",
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/TopTagger.cfg",
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/sampleCollections.cfg",
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/sampleSets.cfg",
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/%(trainingFile)s"%{"trainingFile":mvaFileName},
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/puppiCorr.root",
                      environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so",
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/data/allINone_bTagEff.root", 
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/ISR_Root_Files/ISRWeights.root", 
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/ISR_Root_Files/allINone_ISRJets.root", 
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/CSVv2_Moriond17_B_H.csv", 
                      environ["CMSSW_BASE"] + "/src/TopTaggerTools/Tools/test/data/PileupHistograms_0121_69p2mb_pm4p6.root", 
                      "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_core.so.3.1",
                      "/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/opencv/lib/libopencv_ml.so.3.1",
                      ]

#go make TTopTagger plots!
submitFileTT = """universe = vanilla
Executable = $ENV(CMSSW_BASE)/src/TopTaggerTools/Tools/test/condor/goTTplots.sh
Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $ENV(CMSSW_BASE)/src/TopTaggerTools/Tools/test/condor/goTTplots.sh,TT.tar.gz,$ENV(CMSSW_VERSION).tar.gz
Output = logs/TT_$(Process).stdout
Error = logs/TT_$(Process).stderr
Log = logs/TT_$(Process).log
notify_user = ${LOGNAME}@FNAL.GOV
x509userproxy = $ENV(X509_USER_PROXY)

"""

parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-n',  dest='numfile',               type='int',          default = 5,     help="number of files per job")
parser.add_option ('-d',  dest='datasets',              type='string',       default = '',    help="List of datasets 'ZJetsToNuNu,GJets,DYJetsToLL'")
parser.add_option ('-l',  dest='dataCollections',       action='store_true', default = False, help="List all datacollections")
parser.add_option ('-L',  dest='dataCollectionslong',   action='store_true', default = False, help="List all datacollections and sub collections")
parser.add_option ('-r',  dest='refLumi',               type='string',       default = None,  help="Data collection to define lumi (uses default lumi if no reference data collection is defined)")
parser.add_option ('-c',  dest='noSubmit',              action='store_true', default = False, help="Do not submit jobs.  Only create condor_submit.txt.")

options, args = parser.parse_args()

submitFile = ""
exeName = ""

def makeExeAndFriendsTarrball(filestoTransfer, fname):
    if not options.dataCollections and not options.dataCollectionslong:
        #WORLDSWORSESOLUTIONTOAPROBLEM
        system("mkdir -p WORLDSWORSESOLUTIONTOAPROBLEM")
        for fn in filestoTransfer:
            system("cd WORLDSWORSESOLUTIONTOAPROBLEM; ln -s %s" % fn)
        
        print "Create tarball {0}.tag.gz".format(fname)
        tarallinputs = "tar czvf %s.tar.gz WORLDSWORSESOLUTIONTOAPROBLEM --dereference" % fname
        print tarallinputs
        system(tarallinputs)
        system("rm -r WORLDSWORSESOLUTIONTOAPROBLEM")

submitFile = submitFileTT
fileParts = [submitFile]
sc = SampleCollection("../sampleSets.cfg", "../sampleCollections.cfg")
datasets = []

if options.dataCollections or options.dataCollectionslong:
    scl = sc.sampleCollectionList()
    for sampleCollection in scl:
        sl = sc.sampleList(sampleCollection)
        print sampleCollection
        if options.dataCollectionslong:
            sys.stdout.write("\t")
            for sample in sl:
                sys.stdout.write("%s  "%sample[1])
            print ""
            print ""
    exit(0)

if options.datasets:
    datasets = options.datasets.split(',')
else:
    print "No dataset specified"
    exit(0)

now = datetime.datetime.now()

dirName = "submission_%s"%now.strftime("%Y-%m-%d_%H-%M-%S")
os.system("mkdir %s"%dirName)
os.chdir(dirName)

if not options.dataCollections and not options.dataCollectionslong:
    print "Create tarball ${CMSSW_VERSION}.tar.gz"
    system("tar --exclude-caches-all --exclude-vcs -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp")

# makeExeAndFriendsTarball() is necessary now to apply WORLDSWORSESOLUTIONTOAPROBLEM 

exeName = "simpleAnalyzer"
makeExeAndFriendsTarrball(filestoTransferTT, "TT")

nFilesPerJob = options.numfile

lumis = sc.sampleCollectionLumiList()
lumi = sc.getFixedLumi()
if options.refLumi != None:
    lumi = lumis[options.refLumi]
    print "Normalizing to %s pb-1" % (lumi)

for ds in datasets:
    ds = ds.strip()

    print ds
    # s: file, n:name, e:nEvts
    for s, n, e in sc.sampleList(ds):
        print "\t%s"%n
        #print "\t{0} {1} {2}".format(s, n, e)
        try:
            f = open(s)
        except IOError:
            fShort = s.split("/")[-1]
            if(os.path.isfile(fShort)):
                os.remove(fShort)
            system("xrdcp root://cmseos.fnal.gov/$(echo %s | sed 's|/eos/uscms||') ."%s)
            print "fShort = {0}".format(fShort)
            f = open(fShort)
        if not f == None:
            count = 0
            for l in f:
                if '.root' in l and not 'failed' in l:
                    count = count + 1
            for startFileNum in xrange(0, count, nFilesPerJob):
                fileParts.append("Arguments = %s $ENV(CMSSW_VERSION) %i %i %f %s\n"%(n, nFilesPerJob, startFileNum, lumi, s))
                fileParts.append("Output = logs/%s_%s_%i.stdout\n"%(exeName, n, startFileNum))
                fileParts.append("Error = logs/%s_%s_%i.stderr\n"%(exeName, n, startFileNum))
                fileParts.append("Log = logs/%s_%s_%i.log\n"%(exeName, n, startFileNum))
                fileParts.append("Queue\n\n")

            f.close()

fout = open("condor_submit.txt", "w")
fout.write(''.join(fileParts))
fout.close()

if not options.noSubmit: 
    system('mkdir -p logs')
    system("echo 'condor_submit condor_submit.txt'")
    system('condor_submit condor_submit.txt')

