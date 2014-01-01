# $Id: Makefile,v 1.1 2011/10/17 15:18:39 sigamani Exp $
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

BINFILES =  fitApp.cc

PROGRAMS = $(patsubst %.cc,%,$(BINFILES))


# --- External configuration ---------------------------------
CC         = g++ -Wno-deprecated
CCFLAGS    =  -g 
MFLAGS     = -MM
INCLUDES   =
WORKDIR    = tmp/
LIBDIR     = $(WORKDIR)
OBJDIR     = $(WORKDIR)/objects/
# -------------------------------------------------------------

## su lxplus, T2, pccmsrm*  
ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')
ROOFIT_LIBS := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIB= | sed 's|LIB=||')
ROOFIT_LIBS += $(shell cd $(CMSSW_BASE); scram tool info roofit | grep LIB= | sed 's|LIB=||') 


# su laptop
#ROOFIT_INCLUDE := $(ROOTSYS)/include
#ROOFIT_LIBDIR := $(ROOTSYS)/lib



INCLUDES += -I.  -I.. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/
ROOTSYS  ?= ERROR_RootSysIsNotDefined

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -lRIO

# CC files excluding the binaries
CCFILES=$(filter-out $(BINFILES),$(wildcard *.cc))

# List of all object files to build
OLIST=$(patsubst %.cc,$(OBJDIR)/%.o,$(CCFILES))

# Implicit rule to compile all classes
$(OBJDIR)/%.o : %.cc
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)

$(OBJDIR)/%.o : %.cxx
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)


$(PROGRAMS) : $(OLIST)
	@echo "Linking $@"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS)  $(INCLUDES) $(OLIST) \
	$(ROOTLIBS) $(EXTRALIBS) -o $(WORKDIR)/$@   $(patsubst %,%.cc,$@)

default : analysis 

all : ${PROGRAMS}

clean:
	rm -Rf $(WORKDIR)/*
	@#rm -f $(LIBFILE)
	@rm -Rf *.o

veryclean:
	rm -Rf $(WORKDIR)

#.PHONY : shlib lib default clean
