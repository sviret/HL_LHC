#########################
#
# Configuration file for simple MBias events
# production in tracker only
#
# UE tuning is UE_P6S1
#
# Instruction to run this script are provided on this page:
#
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto
#
# Look at STEP II
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 12/04/2013
#
# Script tested with release CMSSW_6_2_0_SLHC20
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Special geometry (Tracker only)
process.load('DataProduction.SkimGeometry.Sim_SKIM_cff')
process.load('DataProduction.SkimGeometry.GeometryExtendedPhase2TkBEReco_SKIM_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
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


# Generate particle gun events


process.generator = cms.EDFilter("Pythia6GeneratorFilter",
				 pythiaHepMCVerbosity = cms.untracked.bool(False),
				 maxEventsToPrint = cms.untracked.int32(0),
				 pythiaPylistVerbosity = cms.untracked.int32(1),
				 filterEfficiency = cms.untracked.double(1.0),
				 crossSection = cms.untracked.double(79150000000),
				 comEnergy = cms.double(14000.0),
				 PythiaParameters = cms.PSet(
		pythiaUESettings = cms.vstring(
			'MSTU(21)=1 ! Check on possible errors during program execution',
			'MSTJ(22)=2 ! Decay those unstable particles',
			'PARJ(71)=10 . ! for which ctau 10 mm',
			'MSTP(33)=0 ! no K factors in hard cross sections',
			'MSTP(2)=1 ! which order running alphaS',
			'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)',
			'MSTP(52)=2 ! work with LHAPDF',
			'PARP(82)=1.9096 ! pt cutoff for multiparton interactions',
			'PARP(89)=1800. ! sqrts for which PARP82 is set',
			'PARP(90)=0.2479 ! Multiple interactions: rescaling power',
			'MSTP(95)=6 ! CR (color reconnection parameters)',
			'PARP(77)=0.6646 ! CR',
			'PARP(78)=0.5454 ! CR',
			'PARP(80)=0.1 ! Prob. colored parton from BBR',
			'PARP(83)=0.8217 ! Multiple interactions: matter distribution parameter',
			'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter',
			'PARP(62)=1.025 ! ISR cutoff',
			'MSTP(91)=1 ! Gaussian primordial kT',
			'PARP(93)=10.0 ! primordial kT-max',
			'MSTP(81)=21 ! multiple parton interactions 1 is Pythia default',
			'MSTP(82)=4 ! Defines the multi-parton model',
            'PARJ(1) =  0.08 ! HAD diquark suppression',
            'PARJ(2) = 0.21 ! HAD strangeness suppression',
            'PARJ(3) = 0.94 ! HAD strange diquark suppression',
            'PARJ(4) = 0.04 ! HAD vectior diquark suppression',
            'PARJ(11) = 0.35 ! HAD P(vector meson), u and d only',
            'PARJ(12) = 0.35 ! HAD P(vector meson) contains ',
            'PARJ(13) = 0.54 ! HAD P(vector meson), heavy quarks',
            'PARJ(21) = 0.34 ! HAD fragmentation pt',
            'PARJ(25) = 0.63 ! HAD eta0 suppression',
            'PARJ(26) = 0.12 ! HAD eta0 suppression'
			),
	processParameters = cms.vstring('MSEL=0         ! User defined processes',
					'MSUB(11)=1     ! Min bias process',
					'MSUB(12)=1     ! Min bias process',
					'MSUB(13)=1     ! Min bias process',
					'MSUB(28)=1     ! Min bias process',
					'MSUB(53)=1     ! Min bias process',
					'MSUB(68)=1     ! Min bias process',
					'MSUB(92)=1     ! Min bias process, single diffractive',
					'MSUB(93)=1     ! Min bias process, single diffractive',
					'MSUB(94)=1     ! Min bias process, double diffractive',
					'MSUB(95)=1     ! Min bias process'),
         parameterSets = cms.vstring('pythiaUESettings',
				     'processParameters')
	)
)


# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('MBias_100_trOnly_test.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

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
