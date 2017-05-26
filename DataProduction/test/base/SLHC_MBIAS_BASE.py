flat=TILTEDORFLAT

import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

if flat:
	print 'You choose the flat geometry'
	process.load('L1Trigger.TrackTrigger.TkOnlyFlatGeom_cff') # Special config file for TkOnly geometry
else:
	print 'You choose the tilted geometry'
	process.load('L1Trigger.TrackTrigger.TkOnlyTiltedGeom_cff') # Special config file for TkOnly geometry

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEVTS)
)

# Input source
process.source = cms.Source("EmptySource")


# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MYGLOBALTAG', '')


process.RandomNumberGeneratorService.generator.initialSeed      = NSEEDA
process.RandomNumberGeneratorService.VtxSmeared.initialSeed     = NSEEDB
process.RandomNumberGeneratorService.g4SimHits.initialSeed      = NSEEDC

# From
# http://cmslxr.fnal.gov/lxr/source/Configuration/Generator/python/MinBias_14TeV_pythia8_TuneCUETP8M1_cfi.py?v=CMSSW_8_1_0_pre7 
#

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
                          crossSection = cms.untracked.double(71.39e+09),
#                          crossSection = cms.untracked.double(79.15e+09),
                          maxEventsToPrint = cms.untracked.int32(0),
                          pythiaPylistVerbosity = cms.untracked.int32(1),
                          filterEfficiency = cms.untracked.double(1.0),
                          pythiaHepMCVerbosity = cms.untracked.bool(False),
                          comEnergy = cms.double(14000.0),
                          PythiaParameters = cms.PSet(
         pythia8CommonSettingsBlock,
         pythia8CUEP8M1SettingsBlock,
         processParameters = cms.vstring(
             'SoftQCD:nonDiffractive = on',
             'SoftQCD:singleDiffractive = on',
             'SoftQCD:doubleDiffractive = on'),
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
    fileName = cms.untracked.string('MBiasSample.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

#process.g4SimHits.TrackerSD = cms.PSet(
#        ZeroEnergyLoss = cms.bool(False),
#        PrintHits = cms.bool(False),
#        ElectronicSigmaInNanoSeconds = cms.double(12.06),
#        NeverAccumulate = cms.bool(False),
#        EnergyThresholdForPersistencyInGeV = cms.double(0.001),
#        EnergyThresholdForHistoryInGeV = cms.double(0.001)
#    )


# Path and EndPath definitions
process.generation_step      = cms.Path(process.pgen)
process.simulation_step      = cms.Path(process.psim)
process.genfiltersummary_step= cms.EndPath(process.genFilterSummary)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)


process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq

