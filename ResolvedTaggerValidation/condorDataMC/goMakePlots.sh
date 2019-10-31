#!/bin/bash

_PWD=${PWD}

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

tar -zxvf $2.tar.gz
cd $2/
mkdir -p src
cd $2/src
scram b ProjectRename
eval `scramv1 runtime -sh`

#set up local code
tar -xzf ${_CONDOR_SCRATCH_DIR}/DataMC.tar.gz
cd WORLDSWORSESOLUTIONTOAPROBLEM

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}


xrdcp root://cmseos.fnal.gov/$(echo $5 | sed 's|/eos/uscms||') .

ls -lhrt

./DataMC $1 -1 $3 $4 "condor"

ls -lhrt

mv *DataMC*.root ${_CONDOR_SCRATCH_DIR}

rm $(echo $5 | sed 's|.*/||')
rm -r ${_CONDOR_SCRATCH_DIR}/$2