import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("EenTest")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

# process.load("Configuration.StandardSequences.Services_cff")
# #process.load('Configuration.StandardSequences.Geometry_cff')
# #process.load('Configuration/StandardSequences/MagneticField_38T_cff')
# process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
# #process.load("Configuration.StandardSequences.Reconstruction_cff")
# process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
# process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')


process.load('Configuration/EventContent/EventContent_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc'   , '')


#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

#For the non-IC:
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'

#For the IC:
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_IdealEcalIC_v0'

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(900) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 

readFiles.extend([
    # This is an AODSIM example in 80X
    #'file:DoublePhoton_AODSIM_example.root',
    'file:DoubleElectron_AODSIM_example.root',
    ])
secFiles.extend([
    ])

process.source = cms.Source(
    "PoolSource",
    fileNames = readFiles,
    secondaryFileNames = secFiles
    )


########################################
# Define the analyzer
########################################

process.een_analyzer = cms.EDAnalyzer(
    'SimpleNtuplizer',
    vertices     = cms.InputTag("offlinePrimaryVertices"),
    electrons    = cms.InputTag("gedGsfElectrons"),
    photons      = cms.InputTag("photons"),
    rho          = cms.InputTag("fixedGridRhoFastjetAll"),
    genparticles = cms.InputTag("genParticles"),
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('output.root')
    )

process.p = cms.Path(process.een_analyzer)