CXX           = g++
CXXFLAGS      = $(OPT2) -m64 -pipe -Wall -W -Woverloaded-virtual -g
LD            = g++
LDFLAGS       = $(OPT2) -m64 -bind_at_load -L.
SOFLAGS       = -dynamiclib -single_module -undefined $(UNDEFOPT)

ROOTINCLUDES = -I/uscmst1/prod/sw/cms/slc5_amd64_gcc434/lcg/root/5.22.00d-cms30/include
ROOTLDLIBS = -L/uscmst1/prod/sw/cms/slc5_amd64_gcc434/lcg/root/5.22.00d-cms30/lib -lMathCore -lTree -lCore -lNet -lRIO -lCint -lReflex -lCintex -lHist -lGraf -lGraf3d -lGpad  -lRint -lPostscript -lMatrix -lPhysics -lEG -lThread -lm -ldl

INCLUDES = $(PYINCLUDES) $(FASTJETINCLUDES) $(ROOTINCLUDES)
LDLIBS = $(PYLDLIBS) $(FASTJETLDLIBS) $(ROOTLDLIBS)


all : templateCheck
	@echo Build $@

clean:
	rm -rf templateCheck

templateCheck : templateCheck.C
	$(CXX) $(CXXFLAGS) $(INCLUDES)   $^ $(LDLIBS) -o templateCheck


.cc.a :
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^	
	ar rv $@ $*.o
	rm -f $*.o
	@echo Compiled $@
