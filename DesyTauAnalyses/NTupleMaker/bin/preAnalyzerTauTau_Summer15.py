#!/usr/bin/env python

import os
import re
import shlex
import string
import subprocess

#WorkdirLoc = '/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/Sync2015/'
#WorkdirLoc = '/nfs/dust/cms/user/rasp/ntuples/'
WorkdirLoc = '/nfs/dust/cms/group/susy-desy/Run2/MC/Stau/MC_Spring15_25ns_v1/'
OutDir     = '/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/ntuples74/'

options = {
        ###mu+tau samples

        ##DATA
#    'Run2012A-22Jan2013-Data' : {
#    'inputFilePath'  : WorkdirLoc+'Data/Data_2012A_ReReco22Jan_HTT_25July2013_Trees_MuTau_v1/',
#    'outputFileName' : OutDir+'nTupleRun2012A-22Jan2013-Data_MuTau.root',
#    'sample'         : 'Run2012A-22Jan2013-Data',
#    'xSection'       : 0,
#    'skimEff'        : 0,
#    'iJson'          : 7
#    },
#    ##Bkg MC
#    'DYJets_TauTau' : {
#    'inputFilePath'  : WorkdirLoc+'BackgroundsMC/DYJets-50-madgraph-PUS10_MC_Bkg_HTT_25July2013_Trees_MuTau_v4/',
#    'outputFileName' : OutDir+'nTupleDYJets_TauTau_MuTau.root',
#    'sample'         : 'DYJets_TauTau',
#    'xSection'       : 3504,
#    'skimEff'        : 1.0 * 0.317439 * 2474447./9669034,
#    'iJson'          : -1
#    } ,
#    ##Higgs MC
#    'GGFH125' : {
#    'inputFilePath'  : WorkdirLoc+'HiggsSM/GluGluToHToTauTau_M-125_MC_TauTau_v2',
#    'outputFileName' : OutDir+'nTupleGGFH125_TauTau.root',
#    'sample'         : 'GGFH125',
#    'xSection'       : 36.80 ,
#    'skimEff'        : 1.0,
#    'iJson'          : -1,
#    'iDiv'           : 0,                                                      
#    'nDiv'           : 1
#    },
#    'VBFH125' : {
#        'inputFilePath'  : WorkdirLoc+'HiggsSM/VBFHToTauTau_M-125_MC_TauTau_v4',
#        #'inputFilePath'  : WorkdirLoc+'VBF_HToTauTau_M-125_13TeV-powheg-pythia6',
#        'outputFileName' : OutDir+'nTupleVBFH125_TauTau.root',
#        'sample'         : 'VBFH125',
#        'xSection'       : 36.80 ,
#        'skimEff'        : 1.0,
#        'iJson'          : -1,
#        'iDiv'           : 0,
#        'nDiv'           : 1
#        },
    ######MSSM
    'SUSYGGH160' : {
    #'inputFilePath'  : WorkdirLoc+'MC_Spring15_v1/SUSYGluGluToHToTauTau_M-160_PY8_25ns/',
    'inputFilePath'  : WorkdirLoc+'SUSYGluGluToHToTauTau_M-160/',    
    'outputFileName' : OutDir+'nTupleSUSYGGH160_TauTau.root',
    'sample'         : 'SUSYGGH160',
    'xSection'       : 1.0,
    'skimEff'        : 1.0,
    'iJson'          : -1,
    'iDiv'           : 0,
    'nDiv'           : 1
    },
    
}

config_template = string.Template('''
import FWCore.ParameterSet.Config as cms

import os

process = cms.PSet()

process.fwliteInput = cms.PSet(
fileNames = cms.vstring(),
maxEvents = cms.int32(-1),
#outputEvery = cms.uint32(1000)
)

inputFilePath = '$inputFilePath'
inputFiles = os.listdir(inputFilePath)
process.fwliteInput.fileNames = cms.vstring([ os.path.join(inputFilePath, inputFile) for inputFile in inputFiles ])

process.fwliteOutput = cms.PSet(
fileName = cms.string('$outputFileName')
)

process.preAnalyzerTauTau = cms.PSet(
sample = cms.string('$sample'),
analysis = cms.string('$analysis'),
xSection = cms.double($xSection),
skimEff = cms.double($skimEff),
iJson = cms.int32($iJson),
iDiv = cms.int32($iDiv),
nDiv = cms.int32($nDiv)
)
''')

currentDirectory    = os.getcwd()
submissionDirectory = os.path.join(currentDirectory, "Configs")

for sample, option in options.items():
    for analysis in [ 'nominal' ]: #, 'TauUp', 'TauDown' ]:
        if  re.search("Data",sample)!=None and analysis != 'nominal':
            continue
        configOptions = option.copy()
        configOptions['analysis'] = analysis
        if  re.search("Data",sample)==None :
            configOptions['outputFileName'] = configOptions['outputFileName'].replace('.root', '_%s.root' % analysis)
        configFileName = "preAnalyzerTauTau_Summer15_%s_%s_cfg.py" % (sample,analysis)
        configFileName_full = os.path.join(submissionDirectory, configFileName)
        configFile = open(configFileName_full, 'w')
        configConfig = config_template.substitute(configOptions)
        configFile.write(configConfig)
        configFile.close()
