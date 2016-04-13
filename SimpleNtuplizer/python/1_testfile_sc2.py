import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend([
       '/store/mc/RunIIFall15MiniAODv1/DoubleElectron_FlatPt-1To300/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/00BAF127-3EA7-E511-8026-24BE05C6E711.root',
       ])

secFiles.extend([
       ])