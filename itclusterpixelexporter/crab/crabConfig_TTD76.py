from CRABClient.UserUtilities import config
config = config()
dataset_name = "TTBar14TeV_CMSSW_11_3_1_patch1_D76_PU200"

config.section_("General")
config.General.requestName = dataset_name

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ["../data"]
# Name of the CMSSW configuration file
config.JobType.psetName = '../python/ITdigiExporter.py'
config.JobType.pyCfgParams  = ['dataset='+dataset_name]

config.section_("Data")
config.Data.inputDataset = '/TTbar_TuneCP5_14TeV-pythia8/Phase2Spring21DRMiniAOD-PU200Phase2D76_113X_mcRun4_realistic_v7_ext1-v1/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 70
config.Data.totalUnits = -1
config.Data.publication = False
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'ITdigiExporter'
config.Data.outLFNDirBase = '/store/user/syuan/crabjobs/ITdigiExporter/'+config.General.requestName
#config.Data.outLFNDirBase = '/store/user/syuan/crabjobs/ITdigiExporter/20211206'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#    Select input data based on run-ranges

# Where the output files will be transmitted to
config.section_("Site")
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_CH_CERNBOX'
