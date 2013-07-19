#########
#
# Example script to run the python extractor on MC events
# produced with the full geometry
# 
# Usage: cmsRun SLHC_Extr_off.py
#
# More info:
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto
#
# Look at STEP III
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 19/07/2013
#
# Script tested with release CMSSW_6_1_2_SLHC6_patch1
#
#########


import FWCore.ParameterSet.Config as cms

process = cms.Process("MIBextractor")

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi") 
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.Digi_cff')


# Other statements

# Global tag for PromptReco
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# The file you want to extract
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:PGun_example.root'),
                            #fileNames = cms.untracked.vstring('file:PU_10_sample_off.root'),        
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Load the extracto
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

# Tune some options (see MIB_extractor_cfi.py for details)

process.MIBextraction.doMC             = True
process.MIBextraction.doSTUB           = True
process.MIBextraction.doPixel          = True
process.MIBextraction.doL1TT           = True
process.MIBextraction.doMatch          = True

process.MIBextraction.analysisSettings = cms.untracked.vstring(
    "matchedStubs 0",
    "posMatching  1",
    "zMatch  0",
    "maxClusWdth  3",
    "windowSize   -1",
    "pdgSel -1",
    "verbose 0"
    )

process.p = cms.Path(process.MIBextraction)


