import FWCore.ParameterSet.Config as cms

process = cms.Process('EXTR2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("EmptySource")


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Load the extracto
process.load("Extractors.RecoExtractor.MIB_extractor_cff")

process.MIBextraction.fillTree         = False
process.MIBextraction.doBANK           = True
process.MIBextraction.n_events         = 2000

process.MIBextraction.inputRootFile=cms.string('extracted.root')
process.MIBextraction.extractedRootFile=cms.string('extracted_forBK.root')

process.p        = cms.Path(process.MIBextraction)
process.schedule = cms.Schedule(process.p)
