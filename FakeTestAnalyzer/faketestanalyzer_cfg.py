import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

##
## A) Start from a Dataset
##
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#        'file:/bohome/fanfani/Analysis_DiMuons/store/mc/Summer08/Zmumu/GEN-SIM-RECO/IDEAL_V9_v1/0004/24425B53-C988-DD11-BF6A-0015C5E9B2AB.root'
#    )
#)

##
## B) Start from Generator
##
## import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.GlobalTag.globaltag = 'STARTUP_V5::All'
## Generator Input source
process.source = cms.Source("FlatRandomPtGunSource",
    PGunParameters = cms.untracked.PSet(
        MinPhi = cms.untracked.double(-3.14159265359),
        MinPt = cms.untracked.double(0.99),
        PartID = cms.untracked.vint32(-13),
        MaxEta = cms.untracked.double(2.5),
        MaxPhi = cms.untracked.double(3.14159265359),
        MinEta = cms.untracked.double(-2.5),
        MaxPt = cms.untracked.double(1.01)
    ),
    Verbosity = cms.untracked.int32(0),
)
##
##
##

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.demo = cms.EDAnalyzer('FakeTestAnalyzer',
       # output text file  name
       OutTextFileName = cms.untracked.string("OutFileTest.txt"),
       # output text file size (MB), minimal size is 0.35 KB
       OutputSize = cms.untracked.double(5), 
       # sleep time per event (seconds)
       SleepTime = cms.untracked.int32(0)
)

process.p = cms.Path(process.demo)
