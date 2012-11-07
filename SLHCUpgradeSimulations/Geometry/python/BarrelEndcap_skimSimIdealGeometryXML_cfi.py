import FWCore.ParameterSet.Config as cms

### Modifed by Pratima Jindal,Purdue University Calumet, July 2009 to include files for Phase 1 geometry

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('SLHCUpgradeSimulations/Geometry/data/PhaseI/materials.xml', 
        'Geometry/CMSCommonData/data/rotations.xml', 
        'Geometry/CMSCommonData/data/normal/cmsextent.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/cms.xml', 
        'Geometry/CMSCommonData/data/cmsMother.xml', 
        'Geometry/CMSCommonData/data/cmsTracker.xml', 
        'Geometry/CMSCommonData/data/mgnt.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/beampipe.xml', 
        'Geometry/CMSCommonData/data/cmsBeam.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdMaterials.xml', 
        'Geometry/TrackerCommonData/data/pixfwdCommon.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdCylinder.xml', 
        'SLHCUpgradeSimulations/Geometry/data/BarrelEndcap/pixfwd.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdDisks.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdInnerDisk1.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdInnerDisk2.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdInnerDisk3.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdOuterDisk1.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdOuterDisk2.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdOuterDisk3.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdblade1.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdblade2.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixfwdblade3.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarmaterial.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarladder.xml', 
        #'Geometry/TrackerCommonData/data/pixbarladderhalf.xml',
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarladderfull0.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarladderfull1.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarladderfull2.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarladderfull3.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarlayer.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarlayer0.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarlayer1.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarlayer2.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/pixbarlayer3.xml', 
        'SLHCUpgradeSimulations/Geometry/data/BarrelEndcap/pixbar.xml', 
        'Geometry/TrackerCommonData/data/trackermaterial.xml', 
        'SLHCUpgradeSimulations/Geometry/data/BarrelEndcap/tracker.xml', 
        'Geometry/TrackerCommonData/data/trackerpixbar.xml', 
        'SLHCUpgradeSimulations/Geometry/data/PhaseI/trackerpixfwd.xml', 
        'SLHCUpgradeSimulations/Geometry/data/BarrelEndcap/trackerStructureTopology.xml', 
        'SLHCUpgradeSimulations/Geometry/data/BarrelEndcap/trackersens.xml', 
        'SLHCUpgradeSimulations/Geometry/data/BarrelEndcap/trackerRecoMaterial.xml',  
        'SLHCUpgradeSimulations/Geometry/data/BarrelEndcap/trackerProdCuts.xml', 
        'Geometry/TrackerSimData/data/trackerProdCutsBEAM.xml', 
         'Geometry/CMSCommonData/data/FieldParameters.xml'
	
	),
    rootNodeName = cms.string('cms:OCMS')

)


