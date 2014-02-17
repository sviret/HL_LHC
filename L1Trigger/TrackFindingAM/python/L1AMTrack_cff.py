import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.MagneticField_cff import *
from L1Trigger.TrackFindingAM.L1AMTrack_cfi import *
from SimTracker.TrackTriggerAssociation.TTStubAssociation_cfi import * 
from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
import FWCore.ParameterSet.Config as cms


#
# For both steps, in addition to the TTrack container, one produces a container of filtered stubs 
#

#
# STEP 1: AM-based pattern recognition 

TTClusterAssociatorFromPixelDigis.TTClusters  = cms.VInputTag( cms.InputTag("TTPatternsFromStub", "ClusInPattern"))
TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("TTPatternsFromStub", "StubInPattern"))
TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusInPattern"))

TTPatternsFromStubs   = cms.Sequence(TTPatternsFromStub*TTClusterAssociatorFromPixelDigis*TTStubAssociatorFromPixelDigis)


#
# STEP 2: Hough transform fit

#TTClusterAssociatorFromPixelDigis.TTClusters  = cms.VInputTag( cms.InputTag("TTTracksFromPattern", "ClusInTrack"))
#TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("TTTracksFromPattern", "StubInTrack"))
#TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusInTrack"))

TTTracksFromPatterns  = cms.Sequence(TTTracksFromPattern)
