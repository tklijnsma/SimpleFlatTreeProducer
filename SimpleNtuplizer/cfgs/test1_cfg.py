import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("EenTest")

process.load("FWCore.MessageService.MessageLogger_cfi")


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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 

readFiles.extend([
        #'/store/mc/RunIIFall15MiniAODv1/DoubleElectron_FlatPt-1To300/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/00BAF127-3EA7-E511-8026-24BE05C6E711.root',
        #'/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/000B1591-F96E-E411-8885-00266CFFC198.root',
        #'file:/afs/cern.ch/user/i/ikrav/workspace/releases-git/CMSSW_7_2_0/src/ElectronWork/ElectronNtupler/test/00A074A5-BF72-E411-B455-003048F02CBE.root',

        # 80X:

        #'/store/mc/RunIISpring16MiniAODv1/DoublePhoton_FlatPt-5To300/MINIAODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/00000/0231EEB8-9502-E611-BC9C-0025905A6060.root',

        # AODSIM:

        #'/store/mc/RunIISpring16DR80/DoublePhoton_FlatPt-5To300/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/00000/001ED56D-B902-E611-9CCB-000F53273738.root',
        #'/store/mc/RunIISpring16DR80/DoublePhoton_FlatPt-5To300/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/00000/00AFF3C7-7104-E611-A468-02163E00EB4C.root',
        'file:AODSIM_example.root',

       ])
secFiles.extend([
       ])

process.source = cms.Source(
    "PoolSource",
    fileNames = readFiles,
    secondaryFileNames = secFiles
    )


# process.ntupler = cms.EDAnalyzer('SimpleNtuplizer',
#                                  packed = cms.InputTag("packedGenParticles"),
#                                  pruned = cms.InputTag("prunedGenParticles"),
#                                  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                                  pileup = cms.InputTag("addPileupInfo"),
#                                  electrons = cms.InputTag("slimmedElectrons"),
#                                  rho = cms.InputTag("fixedGridRhoFastjetAll")
# )


########################################
# Define the analyzer
########################################

# USE THIS ANALYZER FOR miniAOD SAMPLES

# process.een_analyzer = cms.EDAnalyzer(
#     'SimpleNtuplizer',
#     #'ElectronNtuplerEventStructure',
#     packed    = cms.InputTag("packedGenParticles"),
#     pruned    = cms.InputTag("prunedGenParticles"),
#     vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
#     pileup    = cms.InputTag("addPileupInfo"),
#     electrons = cms.InputTag("slimmedElectrons"),
#     rho       = cms.InputTag("fixedGridRhoFastjetAll")
#     )

# USE THIS ANALYZER FOR AODSIM SAMPLES

process.een_analyzer = cms.EDAnalyzer(
    'SimpleNtuplizer',
    #'ElectronNtuplerEventStructure',
    packed    = cms.InputTag("packedGenParticles"),
    pruned    = cms.InputTag("prunedGenParticles"),
    vertices  = cms.InputTag("offlinePrimaryVertices"),
    pileup    = cms.InputTag("addPileupInfo"),
    electrons = cms.InputTag("gedGsfElectrons"),
    rho       = cms.InputTag("fixedGridRhoFastjetAll")
    )



process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('output.root')
    )

process.p = cms.Path(process.een_analyzer)