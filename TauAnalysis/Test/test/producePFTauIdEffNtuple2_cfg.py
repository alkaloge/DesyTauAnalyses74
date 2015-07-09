import FWCore.ParameterSet.Config as cms

process = cms.Process("producePFTauIdEffNtuple2")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd.ba.infn.it//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
        #'/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/AODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/125A6B71-C56A-E411-9D2B-0025907609BE.root'
        #'file:copyEvent/pickevents_merged.root'
    ),
##    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
##    inputCommands=cms.untracked.vstring(
##        'keep *',
##        'drop patTaus_*_*_*',
##        'drop *PFTau*_*_*_*'
##    )
)
 
#####################################################
  
process.producePFTauIdEffNtuple2Sequence = cms.Sequence()

#--------------------------------------------------------------------------------
# print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(1)
)
process.producePFTauIdEffNtuple2Sequence += process.printGenParticleList
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# rerun tau reconstruction with latest tags

##process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
##process.producePFTauIdEffNtuple2Sequence += process.PFTau
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select "good" reconstructed vertices
process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")

process.producePFTauIdEffNtuple2Sequence += process.selectPrimaryVertex
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# compute event weights for pile-up reweighting
# (Summer'12 MC to 2012 run A data)

##from TauAnalysis.RecoTools.vertexMultiplicityReweight_cfi import vertexMultiplicityReweight
##process.vertexMultiplicityReweight3d2012RunABC = vertexMultiplicityReweight.clone(
##    inputFileName = cms.FileInPath("TauAnalysis/RecoTools/data/expPUpoissonMean_runs190456to208686_Mu17_Mu8.root"),
##    type = cms.string("gen3d"),
##    mcPeriod = cms.string("Summer12_S10")
##)
##process.producePFTauIdEffNtuple2Sequence += process.vertexMultiplicityReweight3d2012RunABC
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")

# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' )
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    jetSource = cms.InputTag('ak5PFJets'),
    jetCorrections = ( 'AK5PF', jetCorrections, "" ),
    outputModules = []
)

# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
   
process.producePFTauIdEffNtuple2Sequence += process.patDefaultSequence
#--------------------------------------------------------------------------------
simpleCutsWP95 = \
    "(userFloat('nHits') <= 999 " + \
    "&& ( (isEB && userFloat('sihih') < 0.01 && userFloat('dPhi') < 0.8 " + \
    "&&            userFloat('dEta') < 0.007 && userFloat('HoE') < 0.15) " + \
    "||   " + \
    "     (isEE && userFloat('sihih') < 0.03 && userFloat('dPhi') < 0.7 " + \
    "&&            userFloat('dEta') < 0.01 && userFloat('HoE') < 999) ))"

process.elecVeto  = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("selectedPatElectrons"),
    cut = cms.string(simpleCutsWP95 + " && pt > 15. && abs(eta) < 2.5 && userFloat('PFRelIsoDB04v3') < 0.3 && abs(userFloat('dzWrtPV')) < 0.2"),
    filter = cms.bool(False)
)
process.producePFTauIdEffNtuple2Sequence += process.elecVeto

process.pfTauIdEffNtuple2Producer = cms.EDProducer("PFTauIdEffNtupleProducer2",
    srcGenParticles = cms.InputTag('genParticles'),
    srcGenJets = cms.InputTag('ak5GenJets'),                                                
    srcRecTaus = cms.InputTag('patTaus'),
    srcRecVetoElectrons = cms.InputTag("elecVeto"),
    srcVertices = cms.InputTag('selectedPrimaryVertexPosition'),
    ##srcWeights = cms.VInputTag('vertexMultiplicityReweight3d2012RunABC')
    srcWeights = cms.VInputTag()                                               
)
process.producePFTauIdEffNtuple2Sequence += process.pfTauIdEffNtuple2Producer

process.p = cms.Path(process.producePFTauIdEffNtuple2Sequence)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("pfTauIdEffNtuple2_Phys14.root")
)

processDumpFile = open('producePFTauIdEffNtuple2.dump', 'w')
print >> processDumpFile, process.dumpPython()




