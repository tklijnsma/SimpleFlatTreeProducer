from CRABClient.UserUtilities import config
config = config()

# 20 June: First test of saturation variables
config.General.requestName              = 'DYToEE_NNPDF30_TreeMaker_v2'

config.General.workArea                 = './'
config.General.transferOutputs          = True
config.General.transferLogs             = True

config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'DY_pset.py'
config.JobType.allowUndistributedCMSSW  = True

config.Data.inputDataset                = '/DYToEE_NNPDF30_13TeV-powheg-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
config.Data.inputDBS                    = 'global'
config.Data.splitting                   = 'FileBased'
config.Data.unitsPerJob                 = 5
config.Data.ignoreLocality              = True
config.Data.publication                 = True

config.Site.storageSite                 = 'T2_US_UCSD'
