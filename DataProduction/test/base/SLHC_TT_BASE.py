flat=TILTEDORFLAT

import FWCore.ParameterSet.Config as cms

process = cms.Process('STUBS')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.Generator.PythiaUESettings_cfi import *
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from L1Trigger.TrackTrigger.TrackTrigger_cff import *

SWTUNING

if flat:
	print 'You choose the flat geometry'
	process.load('L1Trigger.TrackTrigger.TkOnlyFlatGeom_cff') # Special config file for TkOnly geometry
	process.TTStubAlgorithm_official_Phase2TrackerDigi_.zMatchingPS = cms.bool(False) 
	process.TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(
		cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
		cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 2, 3.5, 2, 3.5, 5.5, 6, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 7, 7) ),
		cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 1.5, 3, 2, 3, 5, 6, 6.5, 6.5, 6.5, 5, 6.5, 6.5, 7, 7) ),
		cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 0.5, 0.5, 1, 1.5, 3, 4.5, 6, 6.5, 6.5, 7, 7, 7, 7, 7) ),
		cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 0.5, 0.5, 0.5, 1.5, 2., 3.5, 5., 6.5, 6.5, 6.5, 6, 7, 7, 7) ),
		cms.PSet( EndcapCut = cms.vdouble( 0, 0.5, 0.5, 0.5, 0.5, 1., 1.5, 2.5, 4., 5, 7, 5.5, 7, 7, 7, 7) ),
		)
else:
	print 'You choose the tilted geometry'
	process.load('L1Trigger.TrackTrigger.TkOnlyTiltedGeom_cff') # Special config file for TkOnly geometry

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEVTS)
)

# Input source
process.source = cms.Source("EmptySource")

# The number of pileup events you want  
process.mix.bunchspace=cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.input.nbPileupEvents.averageNumber = cms.double(NPU)             
process.mix.input.fileNames     = cms.untracked.vstring('file:PUFILEA')


# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGLOBALTAG', '')



process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('PGun_example.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)


process.RAWSIMoutput.outputCommands.append('keep  *_*_*_*')
process.RAWSIMoutput.outputCommands.append('drop  *_mix_*_STUBS')
process.RAWSIMoutput.outputCommands.append('drop  PCaloHits_*_*_*')
process.RAWSIMoutput.outputCommands.append('drop  *_ak*_*_*')
process.RAWSIMoutput.outputCommands.append('drop  *_simSi*_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')
process.RAWSIMoutput.outputCommands.append('keep  *_mix_Tracker_*')

process.RandomNumberGeneratorService.generator.initialSeed      = NSEEDA
process.RandomNumberGeneratorService.VtxSmeared.initialSeed     = NSEEDB
process.RandomNumberGeneratorService.g4SimHits.initialSeed      = NSEEDC
process.RandomNumberGeneratorService.mix.initialSeed            = NSEEDD

# for TTbar events
# https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/TTbar_14TeV_TuneCUETP8M1_cfi.py

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         maxEventsToPrint = cms.untracked.int32(0),
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.0),
                         comEnergy = cms.double(14000.0),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
            'Top:gg2ttbar = on ',
            'Top:qqbar2ttbar = on ',
            '6:m0 = 175 ',
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'processParameters',
                                    )
        )
)

# Path and EndPath definitions
process.generation_step         = cms.Path(process.pgen)
process.simulation_step         = cms.Path(process.psim)
process.genfiltersummary_step   = cms.EndPath(process.genFilterSummary)
process.digitisationTkOnly_step = cms.Path(process.pdigi_valid)
process.L1TrackTrigger_step     = cms.Path(process.TrackTriggerClustersStubs)
process.L1TTAssociator_step     = cms.Path(process.TrackTriggerAssociatorClustersStubs)
process.endjob_step             = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step       = cms.EndPath(process.RAWSIMoutput)


process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisationTkOnly_step,process.L1TrackTrigger_step,process.L1TTAssociator_step,process.endjob_step,process.RAWSIMoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq
	

# Automatic addition of the customisation function
from L1Trigger.TrackTrigger.TkOnlyDigi_cff import TkOnlyDigi
process = TkOnlyDigi(process)
# End of customisation functions

