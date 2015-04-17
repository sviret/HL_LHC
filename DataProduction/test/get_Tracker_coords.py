#########
#
# Example script to get the tracker coordinates
# for the skimmed geometry 
# 
# Usage: cmsRun get_Tracker_coords.py
#
# More info:
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto
#
# Look at STEP III
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 16/04/2015
#
# Script tested with release CMSSW_6_2_0_SLHC25_patch3
#
#########


import FWCore.ParameterSet.Config as cms

process = cms.Process("MIBextractor")

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# Special geometry (Tracker only)
process.load('DataProduction.SkimGeometry.GeometryExtendedPhase2TkBEReco_SKIM_cff')

# Other statements

# Global tag for PromptReco
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'DES23_62_V1::All', '')

# Number of events
process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(1)
)

# The file you want to extract

# Input source
process.source = cms.Source("EmptySource")


# Load the extractor
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

# Tune some options (see MIB_extractor_cfi.py for details)

process.MIBextraction.getCoords        = True

process.p = cms.Path(process.MIBextraction)


