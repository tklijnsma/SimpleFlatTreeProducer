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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 

readFiles.extend([
    # This is an AODSIM example in 80X
    # 'file:/nfs-5/users/rclsa/RegressionTraining/TreeMaker/CMSSW_8_0_9/src/SimpleFlatTreeProducer/SimpleNtuplizer/example/SingleElectron_Run2016B_PromptReco_v2.root',
    # 'file:/nfs-5/users/rclsa/RegressionTraining/TreeMaker/CMSSW_8_0_9/src/SimpleFlatTreeProducer/SimpleNtuplizer/example/DYToEE_NNPDF30_13TeV_powheg_pythia8.root'    
    # 'file:/afs/cern.ch/work/r/rcoelhol/public/ExampleAODs/SingleElectron_Run2016B_PromptReco_v2.root',
    # 'file:/afs/cern.ch/work/r/rcoelhol/public/ExampleAODs/DYToEE_NNPDF30_13TeV_powheg_pythia8.root',
    # 'file:/afs/cern.ch/work/t/tklijnsm/EGM/AODexamples/DoubleElectron_AODSIM_example.root',
    'file:/afs/cern.ch/work/t/tklijnsm/EGM/AODexamples/DoublePhoton_AODSIM_example.root',

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
    vertices            = cms.InputTag("offlinePrimaryVertices"),
    electrons           = cms.InputTag("gedGsfElectrons"),
    photons             = cms.InputTag("gedPhotons"),
    rho                 = cms.InputTag("fixedGridRhoFastjetAll"),
    genparticles        = cms.InputTag("genParticles"),
    PUInfoInputTag      = cms.InputTag("addPileupInfo"),
    genEvtInfoInputTag  = cms.InputTag("generator"),

    # caloclusters        = cms.InputTag("caloclusters"),
    # Saturation
    ecalrechitsEB       = cms.InputTag("reducedEcalRecHitsEB"),
    ecalrechitsEE       = cms.InputTag("reducedEcalRecHitsEE"),

    # T&P
    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    HLTTag              = cms.InputTag("TriggerResults","","HLT"),
    HLTObjTag           = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ElecTrig            = cms.untracked.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v*"),
    ElecFilt            = cms.untracked.vstring("hltEle27erWPLooseGsfTrackIsoFilter"),
    
    # isData              = cms.untracked.bool(True)
    isData              = cms.untracked.bool(True)

    )


########################################
# Electron Identification
########################################

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('gedGsfElectrons','','RECO')
setupAllVIDIdsInModule(process,'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',setupVIDElectronSelection)        
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)


# Trigger    
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('output.root')
    )
    
process.p = cms.Path(process.egmGsfElectronIDSequence * process.een_analyzer)
