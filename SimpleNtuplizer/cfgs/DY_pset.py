import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("EenTest")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.load('Configuration/EventContent/EventContent_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc'   , '')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 

readFiles.extend([
    'file:file:000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root',
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
    vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
    electrons           = cms.InputTag("slimmedElectrons"),
    photons             = cms.InputTag("slimmedPhotons"),
    clusters            = cms.InputTag("reducedEgamma:reducedEBEEClusters"),
    rho                 = cms.InputTag("fixedGridRhoFastjetAll"),
    genparticles        = cms.InputTag("prunedGenParticles"),
    PUInfoInputTag      = cms.InputTag("slimmedAddPileupInfo"),
    genEvtInfoInputTag  = cms.InputTag("generator"),

    # caloclusters        = cms.InputTag("caloclusters"),
    # Saturation
    ecalrechitsEB       = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    ecalrechitsEE       = cms.InputTag("reducedEgamma:reducedEERecHits"),

    # T&P
    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    HLTTag              = cms.InputTag("TriggerResults"),
    HLTObjTag           = cms.InputTag("selectedPatTrigger"),
    ElecTrig            = cms.untracked.vstring("HLT_Ele27_WPTight_Gsf_v*"),
    ElecFilt            = cms.untracked.vstring("hltEle27WPTightGsfTrackIsoFilter"),
    
    # isData              = cms.untracked.bool(True)
    isData              = cms.untracked.bool(False)

    )


########################################
# Electron Identification
########################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons',"","PAT")
setupAllVIDIdsInModule(process,'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',setupVIDElectronSelection)        
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)


# Trigger    
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('output.root')
    )
    
process.p = cms.Path(process.egmGsfElectronIDSequence * process.een_analyzer)
