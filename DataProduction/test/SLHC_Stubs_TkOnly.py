#########################
#
# Configuration file for stub production
# production in Tracker-Only geometry 
#
#
# Author: S.Viret (s.viret@ipnl.in2p3.fr)
# Date        : 24/05/2016
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

process = cms.Process('STUBSPROD')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

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
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:PGun_example.root'),
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('EDM_output_with_stubs.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-FEVT')
    )
)

# Additional output definition
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# This line is necessary to keep track of the Tracking Particles

process.RAWSIMoutput.outputCommands.append('keep  *_*_*_*')
process.RAWSIMoutput.outputCommands.append('drop  *_mix_*_STUBS')
process.RAWSIMoutput.outputCommands.append('drop  PCaloHits_*_*_*')
process.RAWSIMoutput.outputCommands.append('drop  *_ak*_*_*')
process.RAWSIMoutput.outputCommands.append('drop  *_simSi*_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')
process.RAWSIMoutput.outputCommands.append('keep  *_mix_Tracker_*')

# Path and EndPath definitions
process.L1TrackTrigger_step   = cms.Path(process.TrackTriggerClustersStubs)
process.L1TTAssociator_step   = cms.Path(process.TrackTriggerAssociatorClustersStubs)
process.endjob_step           = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step     = cms.EndPath(process.RAWSIMoutput)

process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.L1TTAssociator_step,process.endjob_step,process.RAWSIMoutput_step)


