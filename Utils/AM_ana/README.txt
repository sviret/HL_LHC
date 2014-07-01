#############################
#
# AM_ana tool short description
#
# More info on the tutorial page, part 3
#
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto620
#
# Author: S.Viret (viret@in2p3_dot_fr)
# Date  : 01/07/2014
#
#############################

# 
## Introduction
#

The AM_ana tool is part of the /Utils/ package, available in the HL_LHC repository. To get the package, you need to do the following:

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cd $CMSSW_BASE/src   # You need a CMSSW release setup
git-clone git://github.com/sviret/HL_LHC.git
mv HL_LHC/Utils Utils
cd Utils/
tar -zxf tclap.tar.gz
cd AM_ana/
make
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

The AM_ana tool is now compiled, let's see what it can do for you. The code itself is driven by main.cxx and has 3 classes:

--> rates.cxx: stub and cluster rates calculation (part 3.2.2 of the tutorial)
--> efficiencies.cxx: stub and cluster efficiencies calculation (part 3.2.3 of the tutorial)
--> track_eff.cxx: pattern reco and track fit efficiencies (part 6.2.2 of the tutorial)

All the classes cxx/h files are extensively commented. We now provide a rough description of each of them.

One last important note, the AM-ana tool is based on extracted data. To get extracted ROOT files from the EDM ones, you need to use the Reco_Extractor package. Usage of this package is described in part 3.1 of the tutorial.



# 
## Stub rates 
#
   
Command line: ./AM_ana -c rates -i INPUTFILE -o OUTPUTFILE

where INPUTFILE is the directory and name of the extracted ROOT file containing the pile up events to analyze, and OUTPUTFILE the name of the ROOT file in which the rate information will be contained. Rates are computed in Stub/Module/BX (simple counting), for both cluster/stub types (private and official). Rate analysis is then made with the MakeRatesGraphs macro.


# 
## Stub efficiencies
#

Command line: ./AM_ana -c stub_eff -i INPUTFILE -o OUTPUTFILE -t PDGID

where PDGID denote the absolute value of the type of the particle you want to test (the other params are the same than before). The code works as follow. All the MC particles with pT>0.5GeV/c and D0<0.3cm are selected.

If this particles induces a digi in one module, we look if this digi has led to one cluster or one stub in the module. The number of digis is the denominator, and the number of cluster/stubs, the numerator. This is what we call MODULE EFFICIENCY

Then we define the LAYER EFFICIENCY. If there is at least one digi induced by the particle in the layer/disk, we look if there is at leadt one cluster/stub induced by the particle in the layer/disk. In the denominator, we put one whatever the number of layer/disk digis of the particle is, as soon as there is one. We do the same for the stubs/clusters. 

As for the rates, both efficiencies are computed for official/private stubs/clusters

Note: You need to have the tracker digis in, and the TP information


# 
## L1 track reco efficiency
#
 
Command line: ./AM_ana -c L1track_eff -i SECTORFILE -f TESTFILE -o OUTPUTFILE -q D0MAX -p PTMIN -n 1000

where SECTORFILE denotes the name of the csv file obtained from the TkLayout tool. TESTFILE the name of the extracted root file you want to analyze (6.2.1 output), OUTPUTFILE the name of the root file which will contain the result, D0MAX the maximum impact parameter of the tracks you want to select, and PTMIN the minimum pT. n denotes the number of events you want to test, and is automatically truncated to the size of the input file if necessary.

The code is proceeding as follows:

1. It first retrieves the sector definitions, and creates a module to sector table (the csv file contains the opposite).
 
2. Then a loop over the n events is made, and for each events, all the tracking particles satisfying the D0MAX and PTMIN requirements are selected. The denominator is made from all the tracking particles with at least 4 stubs (wherever they are). The other are forgotten forever

3. We then compute in how many sector the particle has at least 4 stubs. This provide SECTOR EFFICIENCY. We also compute the number of distinct layer/disks in which the particle has induced stubs (NHITS)

4. We then loop over the patterns matched for this event, and look how many patterns are containing at least 4 stubs of the particle (in distinct layer/disks). This will give PATTERN EFFICIENCY.

5. We do the same for the fitted tracksn in order to get FIT EFFICIENCY. The global L1 track efficiency if then the product of the 3 others. 

The macro MakePR_EffPlots.C is providing all the summary plots.
 
