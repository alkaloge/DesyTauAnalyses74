import FWCore.ParameterSet.Config as cms

process = cms.Process("copyToCastor")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START70_V7::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    ),
##    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
##    inputCommands = cms.untracked.vstring(
##        'keep *',
##        'drop recoPFTaus_*_*_*'                      
##    )
)

dummyEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('')
    )
)

originalEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_*_*_*'
    )
) 

from Configuration.EventContent.EventContent_cff import *
process.copyToCastorOutputModule = cms.OutputModule("PoolOutputModule",
    AODSIMEventContent,
    fileName = cms.untracked.string(
        'simZplusJets20PUat25ns_13TeV_AOD.root'      
    ),
    maxSize = cms.untracked.int32(1000000000)                                                
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.copyToCastorOutputModule)

