#########################
#
# Configuration file for high pt qcd gun
# production in tracker only
#
# Switch between flat and tilted geometry if provided at the end
#
# Date  : 24/05/2017
#
# Script tested with release CMSSW_10_0_0_pre1
#
#########################
#
# Here you choose if you want flat (True) or tilted (False) geometry
#

flat=False

###################

import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PythiaUESettings_cfi import *

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
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from L1Trigger.TrackTrigger.TrackTrigger_cff import *

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
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("EmptySource")

# Additional output definition

# Global tag for PromptReco
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Random seeds
process.RandomNumberGeneratorService.generator.initialSeed      = 1
process.RandomNumberGeneratorService.VtxSmeared.initialSeed     = 2
process.RandomNumberGeneratorService.g4SimHits.initialSeed      = 3
process.RandomNumberGeneratorService.mix.initialSeed            = 4


# for TTbar events
# https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/TTbar_14TeV_TuneCUETP8M1_cfi.py

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(14000.0),
                         maxEventsToPrint = cms.untracked.int32(0),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
            'HardQCD:all = on',
            'PhaseSpace:pTHatMin = 600.',
            'PhaseSpace:pTHatMax = 800.'
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'processParameters',
                                    )
        )
                         )


# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('QCD_example.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-FEVT')
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


# Path and EndPath definitions
process.generation_step      = cms.Path(process.pgen)
process.simulationTkOnly_step   = cms.Path(process.psim)
process.genfiltersummary_step   = cms.EndPath(process.genFilterSummary)
process.digitisationTkOnly_step = cms.Path(process.pdigi_valid)
process.L1TrackTrigger_step  = cms.Path(process.TrackTriggerClustersStubs)
process.L1TTAssociator_step  = cms.Path(process.TrackTriggerAssociatorClustersStubs)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)


process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulationTkOnly_step,process.digitisationTkOnly_step,process.L1TrackTrigger_step,process.L1TTAssociator_step,process.endjob_step,process.RAWSIMoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq


# Automatic addition of the customisation function
from L1Trigger.TrackTrigger.TkOnlyDigi_cff import TkOnlyDigi
process = TkOnlyDigi(process)
# End of customisation functions
	
