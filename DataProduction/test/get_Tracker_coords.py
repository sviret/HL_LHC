#########
#
# Example script to get the tracker coordinates
# for the skimmed geometry 
# 
# Usage: cmsRun get_Tracker_coords.py
#
# More info:
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTutoDev
#
# Look at STEP III
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 16/04/2015
#
# Script tested with release CMSSW_10_0_0_pre1
#
#########
#
# Here you choose if you want flat (True) or tilted (False) geometry
#

flat=True

###################

import FWCore.ParameterSet.Config as cms

process = cms.Process("MIBextractor")

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')


# Other statements

# Global tag for PromptReco
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

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
process.MIBextraction.flatBarrel       = flat
process.MIBextraction.fullInfo         = False # Do you want all the coords or just the modid

process.p = cms.Path(process.MIBextraction)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms

if flat:
	print 'You choose the flat geometry'
	process.load('L1Trigger.TrackTrigger.TkOnlyFlatGeom_cff') # Special config file for TkOnly geometry
else:
	print 'You choose the tilted geometry'
	process.load('L1Trigger.TrackTrigger.TkOnlyTiltedGeom_cff') # Special config file for TkOnly geometry

# End of customisation functions


