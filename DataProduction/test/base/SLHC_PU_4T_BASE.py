#########################
#
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 7/11/2012
#
#########################

import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PythiaUESettings_cfi import *

process = cms.Process('EXTR')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')

# Special geometry (Tracker only)
process.load('DataProduction.SkimGeometry.Sim_SKIM_cff')
process.load('DataProduction.SkimGeometry.GeometryExtendedPhase2TkBEReco_SKIM_cff')
process.load('DataProduction.SkimGeometry.mixPU_SKIM_cfi')
process.load('DataProduction.SkimGeometry.Digi_SKIM_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEVTS)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/sviret/public/LHE/4tops_SM_10000_events.lhe')
)

# The number of pileup events you want  
process.mix.input.nbPileupEvents.averageNumber = cms.double(NPU)             
process.mix.input.fileNames     = cms.untracked.vstring('file:PUFILEA')


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGLOBALTAG', '')

# Load the extracto
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

process.MIBextraction.doMC             = True
process.MIBextraction.doPixel          = True
process.MIBextraction.doMatch          = True
process.MIBextraction.doSTUB           = True
process.MIBextraction.doL1TT           = True

process.MIBextraction.analysisSettings = cms.untracked.vstring(
    "evtNum RUN",
    "matchedStubs 0",
    "posMatching  1",
    "zMatch  0",
    "maxClusWdth  3",
    "windowSize   -1",
    "thresh THRESHOLD", 
    "pdgSel -1",
    "verbose 0"
    )


process.RandomNumberGeneratorService.generator.initialSeed      = NSEEDA
process.RandomNumberGeneratorService.VtxSmeared.initialSeed     = NSEEDB
process.RandomNumberGeneratorService.g4SimHits.initialSeed      = NSEEDC
process.RandomNumberGeneratorService.mix.initialSeed            = NSEEDD

process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    comEnergy = cms.double(14000.0),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring('MSEL=0         ! User defined processes', 
                        'PMAS(5,1)=4.4   ! b quark mass',
                        'PMAS(6,1)=172.4 ! t quark mass',
			'MSTJ(1)=1       ! Fragmentation/hadronization on or off',
			'MSTP(61)=1      ! Parton showering on or off'),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    ),
    jetMatching = cms.untracked.PSet(
       scheme = cms.string("Madgraph"),
       mode = cms.string("auto"),	# soup, or "inclusive" / "exclusive"
       MEMAIN_etaclmax = cms.double(5.0),
       MEMAIN_qcut = cms.double(30.0),
       MEMAIN_minjets = cms.int32(0),
       MEMAIN_maxjets = cms.int32(12),
       MEMAIN_showerkt = cms.double(0),   # use 1=yes only for pt-ordered showers !
       MEMAIN_nqmatch = cms.int32(5), #PID of the flavor until which the QCD radiation are kept in the matching procedure; 
                                      # if nqmatch=4, then all showered partons from b's are NOT taken into account
				      # Note (JY): I'm not sure what the default is, but -1 results in a throw...
       MEMAIN_excres = cms.string(""),
       outTree_flag = cms.int32(0)        # 1=yes, write out the tree for future sanity check
    )    
)
## so strip simhits are not asked for

process.mergedtruth.simHitCollections = cms.PSet(
        pixel = cms.vstring (
            'g4SimHitsTrackerHitsPixelBarrelLowTof',
            'g4SimHitsTrackerHitsPixelBarrelHighTof',
            'g4SimHitsTrackerHitsPixelEndcapLowTof',
            'g4SimHitsTrackerHitsPixelEndcapHighTof'
        )
    )

# Path and EndPath definitions
process.generation_step      = cms.Path(process.pgen)
process.simulation_step      = cms.Path(process.psim)
process.genfiltersummary_step= cms.EndPath(process.genFilterSummary)
process.digitisation_step    = cms.Path(process.pdigi)
process.L1TrackTrigger_step  = cms.Path(process.L1TrackTrigger)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.p                    = cms.Path(process.MIBextraction)

process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1TrackTrigger_step,process.p)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# Automatic addition of the customisation function
from DataProduction.SkimGeometry.phase2TkCustomsBE_SKIM import customise 

#call to customisation function
process = customise(process)

