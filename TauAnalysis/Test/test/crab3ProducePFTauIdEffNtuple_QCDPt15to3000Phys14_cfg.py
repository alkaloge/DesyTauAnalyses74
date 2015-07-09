from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'FakeRate-Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v2'
config.General.workArea = 'QCD_Pt-15to3000_Tune4C_Flat_13TeV_PU20bx25_trkalmb_PHYS14_25_V1_v1'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'producePFTauIdEffNtuple2_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFN = '/store/user/anayak/TauIdStudies/Phys14/Validation/QCD_11thDec2014/'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
