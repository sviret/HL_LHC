import FWCore.ParameterSet.Config as cms
# Mostly copied from Configuration/StandardSequences/python/Geometry_cff.py
# The R39F16 Version of the Phase 1 Pixel Upgrade
from SLHCUpgradeSimulations.Geometry.BarrelEndcap_skimSimIdealGeometryXML_cfi import *
from Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi import *
# Reconstruction geometry services
#  Tracking Geometry
from Geometry.CommonDetUnit.globalTrackingGeometry_cfi import *

#Tracker
from RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi import *

#  Alignment
from Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometry_cff import *

