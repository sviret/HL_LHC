#######################################
#
# FE analyser package for tracker phase II
#
# Contact: viret_at_in2p3.fr
#
#######################################


The FE analyzer package is providing some tools in order to test the different FE data format and configuration. It uses as input extracted ROOT files containing CMSSW simulated events and produces different output types: text files, rootuples,...

I/ Install and compile

This recipe is for lxplus, if there is a problem with your install, please contact me: viret_at_in2p3_dot_fr

Installation recipe is the following:

[lxplus]> git-clone git://github.com/sviret/FE_ana.git
[lxplus]> cd FE_ana
[lxplus]> make

This creates an executable called FE_ana in the directory, which can for different tasks, described below.

II/ Produce a sequence of trigger data after the CIC chip

./FE_ana -c evtbuild_TRG_CONC  -i INPUT_PHYS -j INPUT_MBIAS -n NBX -f SECTORDEF -o OUTPUT -m BENDCBC -b BENDMPA -k PHYS_PROP -e BLOCK_SIZE

where:

INPUT_PHYS  : extracted root file, or list of extracted root files containing physics events (ie which might pass L1) 
INPUT_MBIAS : extracted root file, or list of extracted root files containing raw min bias events (ie which might not pass L1) 
NBX         :  number of BXs you want to produce
SECTORDEF   : tracker geometry definition (csv file)
OUTPUT      : name of the ouptut ROOT file you want to produce
BENDCBC     : how many bits you want to code stub bend for the CBC (<=5)
BENDMPA     : how many bits you want to code stub bend for the MPA (<=5)
PHYS_PROP   : proportion of physics events in the CIC TRG sequence (in %)
BLOCK_SIZE  : size of the CIC block, in BX units (baseline is 8)

The result of this command is a root file containing, for all the sequences, the content of the TRG block for all the CIC chips. All the words are also stored, in binary formats. This output file can be used as a starting point to produce data patterns (used in parallel with an L1 output)


III/ Produce a sequence of L1 data after the CIC chip

./FE_ana -c evtbuild_L1_CONC  -i INPUT_PHYS -l BLOCK_SIZE -n NBX -r RATE -f SECTORDEF -o OUTPUT

INPUT_PHYS  : extracted root file, or list of extracted root files containing physics events (ie which might pass L1) 
BLOCK_SIZE  : size of the L1 block within the CIC block, in bits (baseline is 64). Just to compute the FIFO size
NBX         : number of BXs you want to produce
RATE        : L1A rate, in kHz. The L1A will be randomly distributed along the BXs, following trigger rules. Average rate will be around RATE
SECTORDEF   : tracker geometry definition (csv file)
OUTPUT      : name of the ouptut ROOT file you want to produce


The result of this command is a root file containing, for all the bunch crossing containing an L1A, the content of the L1 data for all the CIC chips. All the words are also stored, in binary formats. This output file can be used as a starting point to produce data patterns (used in parallel with a TRG output)



IV/ Produce a sequence of trigger data before the CIC chip

./FE_ana -c evtbuild_TRG_FE  -i INPUT_PHYS -j INPUT_MBIAS -n NSEQ -f SECTORDEF -o OUTPUT -k PHYS_PROP

where:

INPUT_PHYS  : extracted root file, or list of extracted root files containing physics events (ie which might pass L1) 
INPUT_MBIAS : extracted root file, or list of extracted root files containing raw min bias events (ie which might not pass L1) 
NSEQ        : number of CIC blocks, of size BLOCK_SIZE (in number of BXs) you want to produce
SECTORDEF   : tracker geometry definition (csv file)
OUTPUT      : name of the ouptut ROOT file you want to produce
PHYS_PROP   : proportion of physics events in the CIC TRG sequence (in %)

The result of this command is a root file containing, for all the sequences, the content of the TRG block for all the MPA/CBC chips (be careful when precising the number of events, these chips are 8 times more than the CIC ones). All the words are also stored, in binary formats. 


V/ Produce a sequence of L1 data before the CIC chip

./FE_ana -c evtbuild_L1_FE   -i INPUT_PHYS -l BLOCK_SIZE -n NBX -r RATE -f SECTORDEF -o OUTPUT

where:

INPUT_PHYS  : extracted root file, or list of extracted root files containing physics events (ie which might pass L1) 
NBX         : number of BXs you want to produce
RATE        : L1A rate, in kHz. The L1A will be randomly distributed along the BXs, following trigger rules. Average rate will be around RATE
SECTORDEF   : tracker geometry definition (csv file)
OUTPUT      : name of the ouptut ROOT file you want to produce

The result of this command is a root file containing, for all the sequences, the content of the TRG block for all the MPA/CBC chips. All the words are also stored, in binary formats.




Command line:
./FE_ana -c pattgen_CONC -i INPUT -o OUTPUT -n NEVT -r RATE

IV/  Collision sequences generation for a complete tracker layer 
