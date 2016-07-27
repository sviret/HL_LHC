import FWCore.ParameterSet.Config as cms

MIBextraction = cms.EDAnalyzer("RecoExtractor",


##
## First you define the name of the  output ROOTfile
##
                               
  extractedRootFile = cms.string('extracted.root'),


##
## Then the name of the input ROOTfile, if you start from already extracted file
##
                               
  inputRootFile     = cms.string('default.root'),
                               

##
## Here you tell is you start from a RECO file (True) or an extracted ROOTuple (False)
##

  fillTree = cms.untracked.bool(True),
                               
## Then the different info you want in (set to false by default)

  # Main stuff                        
  doMC             = cms.untracked.bool(False),          # Extract the MC information (MC tree)
  GenParticles     = cms.InputTag("genParticles", ""),
  TrkParticles     = cms.InputTag("mix" , "MergedTrackTruth"),



  doSTUB           = cms.untracked.bool(False),          # Extract the official STUB information (TkStub tree)
  TTClusters            = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive"),
  TTClustersAssociators = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterInclusive"),
  TTStubs               = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
  TTStubsAssociators    = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),

  CLUS_container   = cms.string( "TTClustersFromPhase2TrackerDigis" ),
  STUB_container   = cms.string( "TTStubsFromPhase2TrackerDigis" ),
  CLUS_name        = cms.string( "ClusterInclusive" ),
  STUB_name        = cms.string( "StubAccepted" ),

  doL1TRK          = cms.untracked.bool(False),          # Extract the official L1track information
#  L1stub_tag       = cms.InputTag( "TTStubsFromPixelDigis", "StubAccepted"), 
  L1pattern_tag    = cms.InputTag( "MergePROutput", "AML1Patterns"),
  L1tc_tag         = cms.InputTag( "MergeTCOutput", "AML1TCs"), 
  L1track_tag      = cms.InputTag( "MergeTCOutput", "AML1TCs"), 
#  L1track_tag      = cms.InputTag( "MergeFITOutput", "AML1Tracks"), 
                               
  # Add Pixel information                              
  doPixel          = cms.untracked.bool(False),          # Extract the Tracker information (Pixel tree)
  digi_tag         = cms.InputTag( "mix","Tracker" ),    # The collection where to find the pixel info
  digisimlink_tag  = cms.InputTag( "mix","Tracker" ),    # The collection where to find the pixel info
  PUinfo_tag       = cms.InputTag( "addPileupInfo","" ), # The collection where to find the PU stuff
  doMatch          = cms.untracked.bool(False),          # Add the simtrack index to each digi                             

                               
  doBANK           = cms.untracked.bool(False),          # Stub container creation for bank (need official stubs)
  doL1TT           = cms.untracked.bool(False),          # Extract the cluster/stub information
  getCoords        = cms.untracked.bool(False),          # Extract the tracker coordinates

  n_events         = cms.untracked.int32(10),            # How many events you want to analyze (only if fillTree=False)
  skip_events      = cms.untracked.int32(0),             # How many events you want to skip (only if fillTree=False)

  # The analysis settings could be whatever you want
  # 
  # Format is "STRING VALUE" where STRING is the name of the cut, and VALUE the value of the cut

  # See demo scripts for usage
  analysisSettings = cms.untracked.vstring()
                              
)
