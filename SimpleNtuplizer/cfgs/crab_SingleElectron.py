from CRABClient.UserUtilities import config
config = config()

# 20 June: First test of saturation variables
config.General.requestName              = 'SingleElectron_2016B_TreeMaker_v2'

config.General.workArea                 = './'
config.General.transferOutputs          = True
config.General.transferLogs             = True

config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'SingleElectron_pset.py'
config.JobType.allowUndistributedCMSSW  = True

config.Data.inputDataset                = '/SingleElectron/Run2016B-PromptReco-v2/AOD'
config.Data.inputDBS                    = 'global'
config.Data.splitting                   = 'LumiBased'
config.Data.unitsPerJob                 = 50
config.Data.lumiMask                    = '/nfs-5/users/rclsa/RegressionTraining/TreeMaker/CMSSW_8_0_9/src/SimpleFlatTreeProducer/SimpleNtuplizer/cfgs/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.ignoreLocality              = True
config.Data.publication                 = True

config.Site.storageSite                 = 'T2_US_UCSD'
