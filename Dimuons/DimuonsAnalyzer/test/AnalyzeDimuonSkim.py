import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process("Dump")

#keep the logging output to a nice level
process.include("FWCore/MessageLogger/data/MessageLogger.cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#input files
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    "file:/bohome/fanfani/CRAB/Zmumu-Summer08-Dimuonskim-v1_1.root"
    )
)
#events to read
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(-1)
    input = cms.untracked.int32(10)
)

#output file name
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("DimuonsPlots.root")
)

process.plotter = cms.EDFilter(
   "DimuonsAnalyzer"
)
## cofiguration to read TnP efficiencies
#process.plotter = cms.EDFilter("DimuonsAnalyzer",
#    #sampletype = cms.untracked.string("Z"),
#    TkEffFile = cms.untracked.string("/bohome/fanfani/Analysis_DiMuons/TnProotfiles/ZW10invpb_muon_TnPType1_eff_test.root"),
#    SAEffFile = cms.untracked.string("/bohome/fanfani/Analysis_DiMuons/TnProotfiles/ZW10invpb_muon_TnPType0_eff_test.root")
#)

process.mypath = cms.Path(process.plotter)

