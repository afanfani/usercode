import FWCore.ParameterSet.Config as cms

process = cms.Process("MUSKIMCSC")


process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/DPGAnalysis/Skims/python/CSCSkim_Tight_cfg.py,v $'),
    annotation = cms.untracked.string('CRAFT CSCSkim Tight')
)


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/mc/Summer09/Zmumu/GEN-SIM-RECO/MC_31X_V3_SD_Mu9-v1/0003/AA9A1B6C-F2AB-DE11-A55B-001EC9D8D48B.root',
        '/store/mc/Summer09/Zmumu/GEN-SIM-RECO/MC_31X_V3_SD_Mu9-v1/0003/7EAEAD53-F0AB-DE11-8F10-00093D121C2E.root'
),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
)

#------------------------------------------
# Load standard sequences.
#------------------------------------------
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V5::All'

process.load("Configuration.StandardSequences.Reconstruction_cff")


#------------------------------------------
# parameters for the CSCSkim module
#------------------------------------------
process.load("CSCSkim_mc_cfi")
process.cscSkim.typeOfSkim = cms.untracked.int32(10)

#### the path

process.CSCSkimTight = cms.Path(process.cscSkim)


#### output 
process.outputSkim = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *','drop *_MEtoEDMConverter_*_*'),
    fileName = cms.untracked.string("CSCSkim_Tight.root"),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('RAW-RECO'),
      filterName = cms.untracked.string('CSCSkim_Tight')
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('CSCSkimTight'))
)

process.outpath = cms.EndPath(process.outputSkim)


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.debugModules.append('CSCSkimTight')
process.MessageLogger.categories.append('CSCSkim')
process.MessageLogger.cout = cms.untracked.PSet(
  threshold     = cms.untracked.string('DEBUG'),
  default       = cms.untracked.PSet( limit = cms.untracked.int32(0) ),
  FwkReport     = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
  CSCSkim    = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
