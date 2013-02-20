import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
    )
# import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.suppressInfo = cms.untracked.vstring( "mkcands" )
process.MessageLogger.suppressWarning = cms.untracked.vstring( "mkcands" )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
MC=False
# Input source
process.source = cms.Source("PoolSource",
                            #skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
   'file:/gpfs_data/local/cms/store/data/Run2012A/MuOnia/AOD/13Jul2012-v1/00000/FCF2DA6C-52CF-E111-99D1-003048678FE0.root'
   ) )
process.source.inputCommands = cms.untracked.vstring(
    "keep *",
    "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__RECO",
    "drop *_MEtoEDMConverter_*_*"
    )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(500) )


#process.load('Configuration.StandardSequences.GeometryExtended_cff') #42x
#process.load("Configuration.StandardSequences.Reconstruction_cff") #42x
process.load('Configuration.Geometry.GeometryIdeal_cff') # 53x

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_P_V41_AN3::All'

process.load('Configuration/EventContent/EventContent_cff')
#
#  Load common sequences
#
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

####################################################################################
##################################good collisions############################################
    
#### 44x
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                                      minimumNDOF = cms.uint32(4) ,
#                                                      maxAbsZ = cms.double(24),
#                                                      maxd0 = cms.double(2)
#                                           )

## 53x                                    
pvSelection = cms.PSet(
  minNdof = cms.double( 4. )
, maxZ    = cms.double( 24. )
, maxRho  = cms.double( 2. )
)

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter", # checks for fake PVs automatically
  filterParams = pvSelection,
  filter       = cms.bool( False ), # use only as producer
  src          = cms.InputTag( 'offlinePrimaryVertices' )
)

process.primaryVertexFilter = process.goodOfflinePrimaryVertices.clone( filter = True ) # filter on good primary vertex


process.noscraping = cms.EDFilter("FilterOutScraping",
applyfilter = cms.untracked.bool(True),
debugOn = cms.untracked.bool(False),
numtrack = cms.untracked.uint32(10),
thresh = cms.untracked.double(0.25)
)


# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)

# Prune generated particles to muons and their parents
process.genMuons = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "drop  *  ",                     # this is the default
            "++keep abs(pdgId) = 13",        # keep muons and their parents
            "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
      )
 )



process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import  addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers
    # with some customization
if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.05
        process.muonMatch.resolveByMatchQuality = True
addDiMuonTriggers(process)
useExistingPATMuons(process,'cleanPatMuons' , addL1Info=False)
changeTriggerProcessName(process, 'HLT')
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
useL1MatchingWindowForSinglets(process)
process.muonL1Info.maxDeltaR     =0.3
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR     = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0
from PhysicsTools.PatAlgos.tools.trackTools import *
makeTrackCandidates(process,                                        #         patAODTrackCands
        label='TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
        tracks=cms.InputTag('generalTracks'), # input track collection
        particleType='pi+',                   # particle type (for assigning a mass)
        preselection='pt > 0.1',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
        selection='pt > 0.1',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
        isolation={},                         # Isolations to use ('source':deltaR; set to {} for None)
        isoDeposits=[],
        mcAs=None                           # Replicate MC match as the one used for Muons
        );                                    # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackCands')
l1cands.addGenMatch = False

process.PATfilter = cms.EDFilter("X3872FilterPAT")


process.mkcands = cms.EDAnalyzer("JPsiPiPiPAT",
     HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
     inputGEN  = cms.untracked.InputTag("genParticles"),
     VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
     DoJPsiMassConstraint = cms.untracked.bool(True),
     SkipPsi2S = cms.untracked.bool(True),
     SameSign = cms.untracked.bool(False),
     DoMonteCarloTree = cms.untracked.bool(False),
     MonteCarloParticleId = cms.untracked.int32(20443),
     MinNumMuPixHits = cms.untracked.int32(1),
     MinNumMuSiHits = cms.untracked.int32(8),
     MaxMuNormChi2 = cms.untracked.double(7),
     MaxMuD0 = cms.untracked.double(10.0),
     MaxJPsiMass = cms.untracked.double(3.25),
     MinJPsiMass = cms.untracked.double(2.95),
     MinNumTrSiHits = cms.untracked.int32(4),
     MinTrPt = cms.untracked.double(0.400),
     JPsiPiPiMaxDR = cms.untracked.double(1.2),
     XCandPiPiMaxDR = cms.untracked.double(1.1),
     UseXDr = cms.untracked.bool(False),
     JPsiPiPiMaxMass = cms.untracked.double(6.0),
     JPsiPiPiMinMass = cms.untracked.double(0.0),
     resolvePileUpAmbiguity = cms.untracked.bool(True),
     addXlessPrimaryVertex = cms.untracked.bool(True),
     Debug_Output = cms.untracked.bool(False),
##
##  use the correct trigger path
##

     TriggersForMatching = cms.untracked.vstring("HLT_DoubleMu4_Jpsi_Displaced_v9","HLT_DoubleMu4_Jpsi_Displaced_v10","HLT_DoubleMu4_Jpsi_Displaced_v11","HLT_DoubleMu4_Jpsi_Displaced_v12"),
     FiltersForMatching = cms.untracked.vstring("hltDisplacedmumuFilterDoubleMu4Jpsi","hltDisplacedmumuFilterDoubleMu4Jpsi","hltDisplacedmumuFilterDoubleMu4Jpsi","hltDisplacedmumuFilterDoubleMu4Jpsi"),
#hltDoubleMu4JpsiDisplacedL3Filtered
     Chi2NDF_Track =  cms.untracked.double(7.0)

)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuOniaRun2012C_Prompt_JPsiPiPiPAT_ntpl.root')
                                   )

# turn off MC matching for the process
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process,['All'],"",None,[])

process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)
## error in 53x, so removing it
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)
process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('onia2MuMuPAT.root'),
        outputCommands = cms.untracked.vstring('drop *',
            #'keep *_genMuons_*_Onia2MuMuPAT',                      # generated muons and parents
            'keep patMuons_patMuonsWithTrigger_*_NTUPLE',    # All PAT muos including general tracks and matches to triggers
        )
)

process.filter = cms.Sequence(process.goodOfflinePrimaryVertices+process.primaryVertexFilter+process.noscraping)
#44x process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)

process.ntup = cms.Path(process.filter*process.patDefaultSequence*process.patMuonsWithTriggerSequence*process.PATfilter*process.mkcands)


process.schedule = cms.Schedule(process.ntup)

