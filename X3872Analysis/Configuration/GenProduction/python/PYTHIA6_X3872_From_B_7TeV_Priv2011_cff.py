import FWCore.ParameterSet.Config as cms

configurationMetadata = cms.untracked.PSet(
                                           version = cms.untracked.string('$Revision: 1.7 $'),
                                           name = cms.untracked.string('$Source: /local/projects/CMSSW/rep/CMSSW/Configuration/GenProduction/python/PYTHIA6_X3872_From_B_7TeV_cff.py,v $'),
                                           annotation = cms.untracked.string('Fall10: Pythia6 generation of non prompt X3872, 7TeV, D6T tune')
                                           )

from Configuration.Generator.PythiaUEZ2Settings_cfi import *

oniafilter = cms.EDFilter("PythiaFilter",
                          Status = cms.untracked.int32(2),   
                          MaxRapidity = cms.untracked.double(1.5),
                          MinRapidity = cms.untracked.double(-1.5),
                          MinPt = cms.untracked.double(7.0),
                          ParticleID = cms.untracked.int32(9120443) ## X(3872) from EvtGen
                          )
generator = cms.EDFilter("Pythia6GeneratorFilter",
                         ExternalDecays = cms.PSet(
                                                   EvtGen = cms.untracked.PSet(
                                                                               operates_on_particles = cms.vint32(0), # 0 (zero) means default list (hardcoded)
                                                                               # you can put here the list of particles (PDG IDs)
                                                                               # that you want decayed by EvtGen
                                                                               use_default_decay = cms.untracked.bool(False),
                                                                               decay_table = cms.FileInPath('GeneratorInterface/ExternalDecays/data/DECAY_NOLONGLIFE.DEC'),
                                                                               particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt.pdl'),
                                                                               user_decay_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/incl_BtoX3872_Jpsipipi.dec'),
                                                                               list_forced_decays = cms.vstring('MyB+','MyB-','MyB0','Myanti-B0'),
                                                                               ),
                                                   parameterSets = cms.vstring('EvtGen')
                                                   ),
                         
                         
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(7000.0),
                         crossSection = cms.untracked.double(48440000000.0),
                         filterEfficiency = cms.untracked.double(0.00013),
                         maxEventsToPrint = cms.untracked.int32(0),
                         PythiaParameters = cms.PSet(
                                                     pythiaUESettingsBlock,
                                                     bbbarSettings = cms.vstring('MSEL = 1'),
                                                     parameterSets = cms.vstring('pythiaUESettings', 
                                                                                 'bbbarSettings')
                                                     )
                         )
mumugenfilter = cms.EDFilter("MCParticlePairFilter",
                             Status = cms.untracked.vint32(1, 1),
                             MaxEta = cms.untracked.vdouble(2.5, 2.5),
                             MinEta = cms.untracked.vdouble(-2.5, -2.5),
                             ParticleCharge = cms.untracked.int32(-1),
                             MinPt = cms.untracked.vdouble(2.5, 2.5),
                             ParticleID1 = cms.untracked.vint32(13),
                             ParticleID2 = cms.untracked.vint32(13)
                             )
ProductionFilterSequence = cms.Sequence(generator*oniafilter*mumugenfilter)



