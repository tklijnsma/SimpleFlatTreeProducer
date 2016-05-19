from CRABClient.UserUtilities import config
config = config()

#config.General.requestName              = 'TKNtup_1805_Photon'
config.General.requestName              = 'TKNtup_1805_Electron'
config.General.workArea                 = './'
config.General.transferOutputs          = True
config.General.transferLogs             = True

config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'test1_cfg.py'
config.JobType.allowUndistributedCMSSW  = True

#config.Data.inputDataset                = '/DoublePhoton_FlatPt-5To300/RunIISpring16DR80-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
config.Data.inputDataset                = '/DoubleElectron_FlatPt-1To300/RunIISpring16DR80-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
config.Data.inputDBS                    = 'global'
config.Data.splitting                   = 'FileBased'
config.Data.unitsPerJob                 = 1
config.Data.ignoreLocality              = True
config.Data.publication                 = True

config.Data.outLFNDirBase               = '/store/user/tklijnsm'

config.Site.storageSite                 = 'T3_CH_PSI'