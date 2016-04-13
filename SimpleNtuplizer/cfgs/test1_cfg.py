import FWCore.ParameterSet.Config as cms

process = cms.Process("EenTest")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 

readFiles.extend([
       '/store/mc/RunIIFall15MiniAODv1/DoubleElectron_FlatPt-1To300/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/00BAF127-3EA7-E511-8026-24BE05C6E711.root',
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

process.een_analyzer = cms.EDAnalyzer(
    'SimpleNtuplizer',
    #'ElectronNtuplerEventStructure',
    packed    = cms.InputTag("packedGenParticles"),
    pruned    = cms.InputTag("prunedGenParticles"),
    vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
    pileup    = cms.InputTag("addPileupInfo"),
    electrons = cms.InputTag("slimmedElectrons"),
    rho       = cms.InputTag("fixedGridRhoFastjetAll")
    )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('output.root')
    )

process.p = cms.Path(process.een_analyzer)