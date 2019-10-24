## Variables set by configure 
OPENCV_DIRECTORY    = 
HDF5_DIRECTORY      = /cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/hdf5/1.8.17-omkpbe2
HDF5_INC_DIRECTORY  = include
TAGGERDIR           = /uscms_data/d3/semrat/CMSSW_9_3_3/src/TopTaggerTools/Tools/test/../../../TopTagger
SUSYANATOOLSDIR     = /uscms_data/d3/semrat/CMSSW_9_3_3/src/TopTaggerTools/Tools/test/../../../SusyAnaTools

CXX       = g++
LD        = g++
CXXFLAGS  = -g -O2 --std=c++11
LDFLAGS   = 
PREFIX    = /usr/local
LDCONFIG  = ldconfig
LIBSUFFIX = so
##End variables set by configure

# directory to put intermediate files 
ODIR       = obj
IDIR       = ../include
SDIR       = ../src
TIDIR      = $(TAGGERDIR)/Tools/cpp
TSDIR      = $(TAGGERDIR)/Tools/cpp
RSDIR      = $(SUSYANATOOLSDIR)/Tools/
RIDIR      = $(SUSYANATOOLSDIR)/Tools

# Top tagger includes
INCLUDESDIRS  += -I. -I$(TAGGERDIR)/..

# Top tagger source diretoriess
TTIDIR     = $(TAGGERDIR)/TopTagger/include
TPIDIR     = $(TAGGERDIR)/CfgParser/include

#openCV include dirs and libs
ifdef OPENCV_DIRECTORY
	INCLUDESDIRS  += -I$(OPENCV_DIRECTORY)/include/ -I$(OPENCV_DIRECTORY)/modules/core/include 
	INCLUDESDIRS  += -I$(OPENCV_DIRECTORY)/modules/video/include -I$(OPENCV_DIRECTORY)/modules/objdetect/include 
	INCLUDESDIRS  += -I$(OPENCV_DIRECTORY)/modules/ml/include/ -I$(OPENCV_DIRECTORY)/modules/photo/include/ 
	INCLUDESDIRS  += -I$(OPENCV_DIRECTORY)/modules/imgproc/include/
	INCLUDESDIRS  += -I$(OPENCV_DIRECTORY)

	LIBS      += -L$(OPENCV_DIRECTORY)/lib -lopencv_ml -lopencv_core
endif

#hdf5 include dir
INCLUDESDIRS  += -I$(HDF5_DIRECTORY)/$(HDF5_INC_DIRECTORY)

## Include ROOT
INCLUDESDIRS  += $(shell root-config --cflags)

## Enable for maximum warning
#CXXFLAGS += -Wall -Wextra -Weffc++

# Flags for generating auto dependancies 
CXXDEPFLAGS = -MMD -MP

#Necessary libraries, the order matters here!
# ROOT libraries
LIBS       = $(shell root-config --libs)

# hdf5 libs
LIBS += -L$(HDF5_DIRECTORY)/lib -lhdf5


# tagger libraries
LIBS      += -L$(TAGGERDIR)/TopTagger/test -lTopTagger

PROGRAMS = simpleAnalyzer fakeRatePlots

all: mkobj $(PROGRAMS)

mkobj:
	@mkdir -p obj

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(IDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(IDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(IDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RIDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RIDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(RSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(RIDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(TSDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(TIDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(TSDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(TIDIR) -I$(INCLUDESDIRS) -o $@ -c $<

$(ODIR)/%.o : $(TSDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -I$(TIDIR) -I$(INCLUDESDIRS) -o $@ -c $<

simpleAnalyzer: $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/TaggerUtility.o $(ODIR)/samples.o $(ODIR)/simpleAnalyzer.o $(ODIR)/customize.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/TTbarCorrector.o $(ODIR)/ISRCorrector.o $(ODIR)/SusyUtility.o
	$(LD) $^ $(LIBS) -o $@

fakeRatePlots: $(ODIR)/fakeRatePlots.o
	$(LD) $^ $(shell root-config --libs) -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d