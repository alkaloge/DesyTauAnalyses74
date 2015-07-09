import FWCore.ParameterSet.Config as cms

process = cms.Process('dumpTriggerInfo')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V15::All')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/user/veelken/CMSSW_5_3_x/skims/data_TauPlusX_2012runD_AOD_1_1_29z.root'                        
#        '/store/user/veelken/CMSSW_5_3_x/skims/simQCDmuEnrichedPt470to600_AOD_1_1_A8V.root'
#       '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
        '/store/user/bluj/TT_Tune4C_13TeV-pythia8-tauola/TTbar_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_100_1_9Gc.root',
        '/store/user/bluj/TT_Tune4C_13TeV-pythia8-tauola/TTbar_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_101_1_dcx.root',
        '/store/user/bluj/TT_Tune4C_13TeV-pythia8-tauola/TTbar_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_102_1_ec9.root',
        '/store/user/bluj/TT_Tune4C_13TeV-pythia8-tauola/TTbar_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_103_1_Qls.root',
        '/store/user/bluj/TT_Tune4C_13TeV-pythia8-tauola/TTbar_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_104_1_Tgc.root',
        '/store/user/bluj/TT_Tune4C_13TeV-pythia8-tauola/TTbar_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_105_1_KeE.root',
        '/store/user/bluj/TT_Tune4C_13TeV-pythia8-tauola/TTbar_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_106_1_tuu.root'

#        'file:/nfs/dust/cms/user/rasp/CMSSW/CMSSW_7_2_2/src/HLT_Mu17_Mu8_SameSign_DPhi_v1.root'

    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.dumpTriggerInfo = cms.EDAnalyzer("DumpTriggerInfo",
    srcL1GtReadoutRecord = cms.InputTag('gtDigis::RECO'),
    srcHLTresults = cms.InputTag('TriggerResults::HLT'),
    srcFilters = cms.InputTag('Tr')                                     
)

process.p = cms.Path(process.dumpTriggerInfo)
