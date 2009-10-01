import FWCore.ParameterSet.Config as cms

#------------------------------------------
# parameters for the CSCSkim module
#------------------------------------------
# For non-cosmic muons, in MC samples
cscSkim = cms.EDFilter(
    "CSCSkim",
    isSimulation       = cms.untracked.bool(True),
    typeOfSkim         = cms.untracked.int32(1),
    rootFileName       = cms.untracked.string('outputDummy.root'),
    histogramFileName  = cms.untracked.string('CSCSkim_histos.root'),
    nLayersWithHitsMinimum  = cms.untracked.int32(3),
    minimumHitChambers      = cms.untracked.int32(3),
    minimumSegments         = cms.untracked.int32(2),
    demandChambersBothSides = cms.untracked.bool(False),
    makeHistograms          = cms.untracked.bool(True),
    whichEndcap = cms.untracked.int32(2),
    whichStation = cms.untracked.int32(3),
    whichRing = cms.untracked.int32(2),
    whichChamber = cms.untracked.int32(24),
#
    cscRecHitTag  = cms.InputTag("csc2DRecHits"),
    cscSegmentTag = cms.InputTag("cscSegments"),
    SAMuonTag     = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    GLBMuonTag    = cms.InputTag("globalMuons"),
    trackTag      = cms.InputTag("generalTracks")
)


