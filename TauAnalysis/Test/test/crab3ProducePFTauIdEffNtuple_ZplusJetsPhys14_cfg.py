from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DYJetsToLL_M-50_13TeV_PU20bx25_PHYS14_25_V1_v3'
config.General.workArea = 'DYJetsLL'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'producePFTauIdEffNtuple2_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
#config.Data.totalUnits = -1
config.Data.outLFN = '/store/user/anayak/TauIdStudies/Phys14/Validation/DYJetsLL_12thDec2014/'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
