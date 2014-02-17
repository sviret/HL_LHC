import FWCore.ParameterSet.Config as cms


TTPatternsFromStub = cms.EDProducer("TrackFindingAMProducer",
   TTInputStubs       = cms.InputTag("TTStubsFromPixelDigis", "StubAccepted"),
   TTInputClusters    = cms.InputTag("TTClustersFromPixelDigis", "ClusterInclusive"),
   TTPatternName      = cms.string("AML1Patterns"),
   TTFiltStubsName    = cms.string("StubInPattern"),
   TTFiltClustersName = cms.string("ClusInPattern"), 
   inputBankFile      = cms.string('/afs/cern.ch/work/s/sviret/testarea/PatternBanks/BE_5D/Eta7_Phi8/ss32_cov40/612_SLHC6_MUBANK_lowmidhig_sec37_ss32_cov40.pbk'),
   threshold          = cms.int32(5)
)

TTTracksFromPattern = cms.EDProducer("TrackFitHoughProducer",
   TTInputStubs       = cms.InputTag("TTPatternsFromStub", "StubInPattern"),
   TTInputClusters    = cms.InputTag("TTPatternsFromStub", "ClusInPattern"),
   TTInputPatterns    = cms.InputTag("TTPatternsFromStub", "AML1Patterns"),
   TTL1TracksName     = cms.string("AML1Tracks"),
   TTL1StubsName      = cms.string("StubInTrack"),
   TTL1ClustersName   = cms.string("ClusInTrack")                                     
)
