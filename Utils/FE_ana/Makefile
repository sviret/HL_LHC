CXX = g++
LD  = g++
CFLAGS = -Wall -g -std=c++11

INCS = $(ROOTSYS)/include/ $(PWD)/../tclap-1.2.1/include/ .

%.o: %.cxx 
	@echo "*"
	@echo "* compile "$@
	@echo "*"
	$(CXX) $(CFLAGS) $(addprefix -I, $(INCS)) -c $< -o $@

FE_ana:main.o evtbuilder.o patterngen.o jobparams.o l1_builder.o stim_builder.o compa_files.o
	@echo "Build FE data analyzer tool" 
	$(LD) $^ $(shell $(ROOTSYS)/bin/root-config --libs) -o $@

all : FE_ana

clean: 
	\rm *.o	
