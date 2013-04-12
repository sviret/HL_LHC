#########
#
# Example script to run the python extractor on MC events
# 
# Usage: cmsRun SLHC_Extr_test.py
#
# More info:
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto
#
# Look at STEP III
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 12/04/2013
#
# Script tested with release CMSSW_6_1_2_SLHC1
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

process.load('DataProduction.SkimGeometry.GeometryExtendedPhase2TkBEReco_SKIM_cff')
process.load('DataProduction.SkimGeometry.Digi_SKIM_cff')

# Other statements

# Global tag for PromptReco
process.GlobalTag.globaltag = 'POSTLS161_V15::All'

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
                            #fileNames = cms.untracked.vstring('file:PU_10_sample.root'),                           
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# Load the extracto
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

# Tune some options (see MIB_extractor_cfi.py for details)

process.MIBextraction.doMC             = True
process.MIBextraction.doPixel          = True
process.MIBextraction.doL1TT           = True
process.MIBextraction.doMatch          = True

process.MIBextraction.analysisSettings = cms.untracked.vstring(
    "matchedStubs 0",
    "posMatching  1",
    "maxClusWdth  3",
    "windowSize   15",
    "pdgSel -1",
    "verbose 0"
    )

process.p = cms.Path(process.MIBextraction)


