#########################
#
# Configuration file for simple PGUN events
# production in tracker only
#
# Instruction to run this script are provided on this page:
#
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto
#
# Look at STEP II
#
# Author: S.Viret (vioret@in2p3.fr)
# Date  : 7/11/2012
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Special geometry (Tracker only)
process.load('SLHCUpgradeSimulations.Geometry.Digi_skimBarrelEndcap_cff')
process.load('SLHCUpgradeSimulations.Geometry.fakeConditions_BarrelEndcap_cff')
process.load('SLHCUpgradeSimulations.Geometry.recoFromSimDigis_BarrelEndcap_cff')
process.load('SLHCUpgradeSimulations.Geometry.BarrelEndcap_skimSimIdealGeometryXML_cff')
process.load('SimG4Core.Application.g4SimHits_JustTrack_cfi')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")


# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'DESIGN42_V11::All'

# Random seeds
process.RandomNumberGeneratorService.generator.initialSeed      = 1
process.RandomNumberGeneratorService.VtxSmeared.initialSeed     = 2
process.RandomNumberGeneratorService.g4SimHits.initialSeed      = 3
process.RandomNumberGeneratorService.mix.initialSeed            = 4

# Generate particle gun events
process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MinPt  = cms.double(10.),
        MaxPt  = cms.double(50.),
        PartID = cms.vint32(-13),
        MinEta = cms.double(-2.5),
        MaxEta = cms.double(2.5),
        MinPhi = cms.double(0.),
	MaxPhi = cms.double(6.28)
    ),
    Verbosity = cms.untracked.int32(0),
    AddAntiParticle = cms.bool(False),
)


# Output definition

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

process.RAWSIMoutput.outputCommands.append('keep PixelDigiedmDetSetVector_simSiPixelDigis_*_*')
process.RAWSIMoutput.outputCommands.append('keep *_simSiPixelDigis_*_*')

process.mergedtruth.simHitCollections.pixel = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                         'g4SimHitsTrackerHitsPixelBarrelHighTof',
                         'g4SimHitsTrackerHitsPixelEndcapLowTof',
                         'g4SimHitsTrackerHitsPixelEndcapHighTof')

# Path and EndPath definitions
process.generation_step      = cms.Path(process.pgen)
process.psim                 = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")*process.g4SimHits)
process.simulation_step      = cms.Path(process.psim)
process.genfiltersummary_step= cms.EndPath(process.genFilterSummary)
process.digitisation_step    = cms.Path(process.pdigi)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.endjob_step,process.RAWSIMoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
