import FWCore.ParameterSet.Config as cms

process = cms.Process('EXTR')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')

# Special geometry (Tracker only)
process.load('DataProduction.SkimGeometry.Sim_SKIM_cff')
process.load('DataProduction.SkimGeometry.GeometryExtendedPhase2TkBEReco_SKIM_cff')
process.load('DataProduction.SkimGeometry.mixNoPU_SKIM_cfi')
process.load('DataProduction.SkimGeometry.Digi_SKIM_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(NEVTS)
)

# Input source
process.source = cms.Source("EmptySource")


# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MYGLOBALTAG'

# Load the extracto
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

process.MIBextraction.doMC             = True
process.MIBextraction.doPixel          = True
process.MIBextraction.doMatch          = True
process.MIBextraction.doSTUB           = True

process.RandomNumberGeneratorService.generator.initialSeed      = NSEEDA
process.RandomNumberGeneratorService.VtxSmeared.initialSeed     = NSEEDB
process.RandomNumberGeneratorService.g4SimHits.initialSeed      = NSEEDC
process.RandomNumberGeneratorService.mix.initialSeed            = NSEEDD


process.generator = cms.EDProducer("FlatRandomOneOverPtGunProducer",
    PGunParameters = cms.PSet(
	MaxOneOverPt = cms.double(PTMAX),
        MinOneOverPt = cms.double(PTMIN),
	XFlatSpread = cms.double(1.5),  # In mm
	YFlatSpread = cms.double(1.5),  # In mm
	ZFlatSpread = cms.double(150.),  # In mm
        PartID = cms.vint32(PTYPE),
        MaxEta = cms.double(ETAMAX),
	MaxPhi = cms.double(PHIMAX),
        MinEta = cms.double(ETAMIN),
        MinPhi = cms.double(PHIMIN)
    ),
    Verbosity = cms.untracked.int32(0),
    AddAntiParticle = cms.bool(False),
)

process.MIBextraction.doL1TT           = False


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
