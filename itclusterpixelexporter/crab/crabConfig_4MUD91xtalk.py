from CRABClient.UserUtilities import config
config = config()
dataset_name = "RelValFourMuExtendedPt1_200_D91_xtalk"

config.section_("General")
config.General.requestName = dataset_name + "_V1"

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ["../data"]
# Name of the CMSSW configuration file
config.JobType.psetName = '../python/ITdigiExporter.py'
config.JobType.pyCfgParams  = ['dataset='+dataset_name]

config.section_("Data")
config.Data.inputDataset = '/RelValFourMuExtendedPt1_200/CMSSW_12_3_0_pre6-PU_123X_mcRun4_realistic_v8_JIRA146_D91_XTalk-v2/GEN-SIM-DIGI-RAW'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = 50
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
