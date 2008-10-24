// -*- C++ -*-
//
// Package:    DimuonsAnalyzer
// Class:      DimuonsAnalyzer
// 
/**\class DimuonsAnalyzer DimuonsAnalyzer.cc Dimuons/DimuonsAnalyzer/src/DimuonsAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/1
//         Created:  Fri May 30 15:05:01 CEST 2008
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>
                                                                                                                               
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
                                                                                                                               
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
                                                                                                                               
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
                                                                                                                               
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//match?
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
// had to pick it up from CMSSW2012
#include "PhysicsTools/HepMCCandAlgos/interface/MCCandMatcher.h"
//#include "DataFormats/Candidate/interface/CompositeRefBaseCandidate.h"

// add Recocand
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"

// add PAT
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/View.h"

// TFileService files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
                                                                                                                              
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObject.h>
#include <TTree.h>
#include <vector>                                                                                                                               
                                                                                                                               
using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;


//
// class decleration
//

class DimuonsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DimuonsAnalyzer(const edm::ParameterSet&);
      edm::InputTag mymatch_,muonCands_;

      ~DimuonsAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      double readefficiency2D(double pt, double eta);
      double readefficiency(double pt, double eta);
      double readTKefficiency(double pt, double eta);
      Particle::LorentzVector getParticleP4(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2); 
//      bool isZmumu(const Candidate *genCand );

      // ----------member data ---------------------------
  int nEvt;// used to count the number of events
       
  bool TwoMuNonIso;
   
  TFile *theefffile;
  TFile *theSAefffile;                                                                                                                      

// gen inv mass
  TH1F *m_massGen;
  TH1F *m_mumu_Gen;
  TH1F *m_massGenZ;   
  TH1F *m_mumu_GenZ;
  TH1F *CandpdgId;

// inv mass
  TH1F *m_mumu;
  TH1F *m_mutrk;
  TH1F *m_musamu;

// inv mass for 2 muons types
  TH1F *m_mumu_os;
  TH1F *m_mumu_ss;
  TH1F *m_mumu_2Global_ss;
  TH1F *m_mumu_GlobalSA_ss;
  TH1F *m_mumu_2SA_ss;
  TH1F *m_mumu_GlobalTkMu_ss;
  TH1F *m_mumu_2TkMu_ss;
  TH1F *m_mumu_GlobalTkMu_trk_ss;
  TH1F *m_mumu_2TkMu_trk_ss;

  TH1F *m_mumu_GlobalCaloMu_ss;
  TH1F *m_mumu_2Global;
  TH1F *m_mumu_GlobalSA;
  TH1F *m_mumu_2SA;
  TH1F *m_mumu_GlobalTkMu;
  TH1F *m_mumu_2TkMu;
  TH1F *m_mumu_GlobalTkMu_trk;
  TH1F *m_mumu_2TkMu_trk;
  TH1F *m_mumu_GlobalCaloMu;
 
  TH1F *m_mumu_GlobalTkMuAND2TkMu;
  TH1F *m_mumu_GlobalTkMuAND2Glb;
  TH1F *m_mumu_GlobalTkMuAND2TkMu_ss;
  TH1F *m_mumu_GlobalTkMuAND2Glb_ss;
  TH1F *m_mumu_GlobalTkSub;
  TH1F *m_mumu_GlobalTkSub_ss;
  TH1F *m_mumu_GlobalTk_TkMuNotGlb;
  TH1F *m_mumu_GlobalTk_TkMuNotGlb_ss;
  TH1F *m_mumu_GlobalTk_GlbNotTkMu;
  TH1F *m_mumu_GlobalTk_GlbNotTkMu_ss;

  TH1F *m_mutrk_GlobalTk_ss;
  TH1F *m_mutrk_GlobalTk;                                                                                                                             
  TH1F *m_mumu_effcor;
  TH1F *m_mumu_os_effcor;
  TH1F *m_mumu_ss_effcor;
  TH1F *m_mumu_2Global_ss_effcor;
  TH1F *m_mumu_2Global_effcor;
  TH1F *m_mutrk_GlobalTk_ss_effcor;
  TH1F *m_mutrk_GlobalTk_effcor;

// eta vs pt
  TH2F *etavspt_mu_Gen;
  TH2F *etavspt_mu;
  TH2F *etavspt_mu_NoMCMatch;
  TH2F *etavspt_mu_Glb;
  TH2F *etavspt_mu_TkMu;
  TH2F *etavspt_mu_SA;
  TH2F *etavspt_patmu;
  TH2F *etavspt_patmu_Glb;
  TH2F *etavspt_patmu_TkMu;
  TH2F *etavspt_patmu_SA;

  TH2F *etavspt_GlobalTkMu_lowm;
  TH2F *etavspt_GlobalTkMu_lowm_ss;
  TH2F *etavspt_GlobalTkMuSub;
  TH2F *etavspt_GlobalTkMuSub_ss;
  TH2F *etavspt_GlobalTk_GlbNotTkMu;
  TH2F *etavspt_GlobalTk_GlbNotTkMu_ss;
  TH2F *etavspt_GlobalTk_TkMuNotGlb;
  TH2F *etavspt_GlobalTk_TkMuNotGlb_ss;

// eta vs pt for mu not matching MC
  TH2F *etavspt_muNoMC;
  TH2F *etavspt_patmuNoMC;
  TH2F *etavspt_patmuNoMC_Glb;
  TH2F *etavspt_patmuNoMC_TkMu;
  TH2F *etavspt_patmuNoMC_SA;

// charge * muon pt
  TH1F *qpt_mu_Gen;
  TH1F *qpt_patmu;
  TH1F *qpt_patmu_Glb;
  TH1F *qpt_patmu_TkMu;
  TH1F *qpt_patmu_SA;
  TH1F *qpt_mu;
  TH1F *qpt_mu_Glb;
  TH1F *qpt_mu_TkMu;
  TH1F *qpt_mu_SA;


// muon pt
  TH1F *pt_mu_Gen;

  TH1F *pt_patmu;
  TH1F *pt_patmu_Glb;
  TH1F *pt_patmu_Glb_TkMu;
  TH1F *pt_patmu_GlbNotTkMu;
  TH1F *pt_patmu_Glb_SA;
  TH1F *pt_patmu_TkMu;
  TH1F *pt_patmu_TkMu_Glb;
  TH1F *pt_patmu_TkMuNotGlb;
  TH1F *pt_patmu_TkMu_SA;
  TH1F *pt_patmu_SA;
  TH1F *pt_patmu_SA_Glb;
  TH1F *pt_patmu_SA_TkMu;
  TH1F *ptgen_patmuMC_Glb;       // pt gen for matched mu 
  TH1F *ptgen_patmuMC_Glb_TkMu;
  TH1F *ptgen_patmuMC_GlbNotTkMu;
  TH1F *ptgen_patmuMC_TkMu;
  TH1F *ptgen_patmuMC_TkMuNotGlb;
  TH1F *ptgen_patmuMC_SA;
  TH1F *ptgen_patmuMC_Glb_SA;
  TH1F *ptgen_patmuMC_TkMu_SA;

  TH1F *Deltapt_patmuMC; // pt PAT reco - pt gen matched
  TH1F *Deltapt_patmuMC_Glb;
  TH1F *Deltapt_patmuMC_Glb_TkMu;
  TH1F *Deltapt_patmuMC_GlbNotTkMu;
  TH1F *Deltapt_patmuMC_TkMu;
  TH1F *Deltapt_patmuMC_TkMuNotGlb;
  TH1F *Deltapt_patmuMC_SA;
  TH1F *Deltapt_patmuMC_TkMu_SA;
  TH1F *Deltapt_patmuMC_Glb_SA;

  TH1F *Deltapt_muMC; // pt reco - pt gen matched
  TH1F *pt_mu;
  TH1F *pt_mu_Glb;
  TH1F *pt_mu_Glb_TkMu;
  TH1F *pt_mu_Glb_SA;
  TH1F *pt_mu_TkMu;
  TH1F *pt_mu_TkMu_Glb;
  TH1F *pt_mu_TkMu_SA;
  TH1F *pt_mu_SA;
  TH1F *pt_mu_SA_Glb;
  TH1F *pt_mu_SA_TkMu;
  
  TH1F *pt_mu_effcor;
  TH1F *pt_low_mu;
  TH1F *pt_low_mu_effcor;
  TH1F *eta_mu;
  TH1F *eta_mu_effcor;
  TH1F *pt_mu_mutrk;
  TH1F *eta_mu_mutrk;

// 2MuNonIso
  TH1F *m_mumu_2mni;
  TH1F *m_mutrk_2mni;
  TH1F *m_musamu_2mni;

  TH1F *m_mumu_effcor_2mni;

  TH1F *pt_mu_2mni;
  TH1F *pt_mu_effcor_2mni;
  TH1F *pt_low_mu_2mni;
  TH1F *pt_low_mu_effcor_2mni;
  TH1F *eta_mu_2mni;
  TH1F *eta_mu_effcor_2mni;
  TH1F *pt_mu_mutrk_2mni;
  TH1F *eta_mu_mutrk_2mni;

                         
// trigger bits                                                                                                      
  TH1F *triggerbits;
  edm::TriggerNames trigNames ;
                                                                                                                               
  int nDimuonCand, nTrkDimuonCand, nSADimuonCand, nMatchedMu;
  int nTrig;
  int DIMUONMAX;// used to set maximum of arrays
  int nPATMuon, nPATMuonMCMatch, MUONMAX;
  double DimuonCand_goodmumu_mass[10];
  double DimuonCand_mumuonetrack_mass[10];
  double DimuonCand_mumuonesamuon_mass[10];
  double MuonCand_pt1[10], MuonCand_pt2[10];
  double MuonCand_eta1[10], MuonCand_eta2[10];
  double Dimuon_SAG_mass;

  int nbinspt;
  double ptlow;
  double pthigh;
  double binwidth;

  string TKFileName_  ;
  string SAFileName_  ;
  string sampletype_ ;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DimuonsAnalyzer::DimuonsAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  nEvt=0;
  DIMUONMAX=10;
  MUONMAX=10;
  nTrig = 4;
                                                                                                                               
  nbinspt = 20; ptlow = 0.0; pthigh = 100.0; binwidth = (pthigh-ptlow)/nbinspt;
 
            
  //TFileService
  edm::Service<TFileService> fs;
   
  TKFileName_ = iConfig.getUntrackedParameter<string> ("TkEffFile", "");
  SAFileName_ = iConfig.getUntrackedParameter<string> ("SAEffFile", "");
  sampletype_ = iConfig.getUntrackedParameter<string> ("sampletype","whatever");   
  
  if ( ! SAFileName_.empty() ) {                                                                          
   theefffile = TFile::Open(TKFileName_.c_str());
  }
  if ( ! TKFileName_.empty() ) {
   theSAefffile = TFile::Open(SAFileName_.c_str());
  }                               
                                                                           
  //theefffile = TFile::Open("/bohome/fanfani/CRAB/muon_eff_test.root");
                                                                                                                               
  m_mumu = fs->make<TH1F>("m_mumu","m_mumu",500,0,200);
  m_mutrk = fs->make<TH1F>("m_mutrk","m_mutrk",500,0,200);
  m_musamu = fs->make<TH1F>("m_musamu","m_musamu",500,0,200);
  
  m_massGen = fs->make<TH1F>("m_massGen","m_massGen",500,0,200);
  m_mumu_Gen = fs->make<TH1F>("m_mumu_Gen","m_mumu_Gen",500,0,200);
  m_massGenZ = fs->make<TH1F>("m_massGenZ","m_massGenZ",500,0,200);
  m_mumu_GenZ = fs->make<TH1F>("m_mumu_GenZ","m_mumu_GenZ",500,0,200);
 
  CandpdgId  = fs->make<TH1F>("CandpdgId","CandpdgId",40000,0,40000);

  m_mumu_os = fs->make<TH1F>("m_mumu_os","m_mumu_os",500,0,200);
  m_mumu_ss = fs->make<TH1F>("m_mumu_ss","m_mumu_ss",500,0,200);
                                                                                                   
  m_mumu_2Global_ss = fs->make<TH1F>("m_mumu_2Global_ss","m_mumu_2Global_ss",500,0,200);
  m_mumu_GlobalSA_ss = fs->make<TH1F>("m_mumu_GlobalSA_ss","m_mumu_GlobalSA_ss",500,0,200);
  m_mumu_2SA_ss = fs->make<TH1F>("m_mumu_2SA_ss","m_mumu_2SA_ss",500,0,200);
  m_mumu_GlobalTkMu_ss = fs->make<TH1F>("m_mumu_GlobalTkMu_ss","m_mumu_GlobalTkMu_ss",500,0,200);
  m_mumu_2TkMu_ss = fs->make<TH1F>("m_mumu_2TkMu_ss","m_mumu_2TkMu_ss",500,0,200);
  m_mumu_GlobalTkMu_trk_ss = fs->make<TH1F>("m_mumu_GlobalTkMu_trk_ss","m_mumu_GlobalTkMu_trk_ss",500,0,200);
  m_mumu_2TkMu_trk_ss = fs->make<TH1F>("m_mumu_2TkMu_trk_ss","m_mumu_2TkMu_trk_ss",500,0,200);
  m_mumu_GlobalCaloMu_ss = fs->make<TH1F>("m_mumu_GlobalCaloMu_ss","m_mumu_GlobalCaloMu_ss",500,0,200);
  m_mumu_2Global = fs->make<TH1F>("m_mumu_2Global","m_mumu_2Global",500,0,200);
  m_mumu_GlobalSA = fs->make<TH1F>("m_mumu_GlobalSA","m_mumu_GlobalSA",500,0,200);
  m_mumu_2SA = fs->make<TH1F>("m_mumu_2SA","m_mumu_2SA",500,0,200);
  m_mumu_GlobalTkMu = fs->make<TH1F>("m_mumu_GlobalTkMu","m_mumu_GlobalTkMu",500,0,200);
  m_mumu_2TkMu = fs->make<TH1F>("m_mumu_2TkMu","m_mumu_2TkMu",500,0,200);
  m_mumu_GlobalTkMu_trk = fs->make<TH1F>("m_mumu_GlobalTkMu_trk","m_mumu_GlobalTkMu_trk",500,0,200);
  m_mumu_2TkMu_trk = fs->make<TH1F>("m_mumu_2TkMu_trk","m_mumu_2TkMu_trk",500,0,200);
  m_mumu_GlobalCaloMu = fs->make<TH1F>("m_mumu_GlobalCaloMu","m_mumu_GlobalCaloMu",500,0,200);

  m_mumu_GlobalTkMuAND2TkMu = fs->make<TH1F>("m_mumu_GlobalTkMuAND2TkMu","m_mumu_GlobalTkMuAND2TkMu",500,0,200);
  m_mumu_GlobalTkMuAND2TkMu_ss = fs->make<TH1F>("m_mumu_GlobalTkMuAND2TkMu_ss","m_mumu_GlobalTkMuAND2TkMu_ss",500,0,200);
  m_mumu_GlobalTkMuAND2Glb = fs->make<TH1F>("m_mumu_GlobalTkMuAND2Glb","m_mumu_GlobalTkMuAND2Glb",500,0,200);
  m_mumu_GlobalTkMuAND2Glb_ss = fs->make<TH1F>("m_mumu_GlobalTkMuAND2Glb_ss","m_mumu_GlobalTkMuAND2Glb_ss",500,0,200);
  m_mumu_GlobalTkSub  = fs->make<TH1F>("m_mumu_GlobalTkSub","m_mumu_GlobalTkSub",500,0,200);
  m_mumu_GlobalTkSub_ss = fs->make<TH1F>("m_mumu_GlobalTkSub_ss","m_mumu_GlobalTkSub_ss",500,0,200);
  m_mumu_GlobalTk_TkMuNotGlb = fs->make<TH1F>("m_mumu_GlobalTk_TkMuNotGlb","m_mumu_GlobalTk_TkMuNotGlb",500,0,200);
  m_mumu_GlobalTk_TkMuNotGlb_ss = fs->make<TH1F>("m_mumu_GlobalTk_TkMuNotGlb_ss","m_mumu_GlobalTk_TkMuNotGlb_ss",500,0,200);
  m_mumu_GlobalTk_GlbNotTkMu = fs->make<TH1F>("m_mumu_GlobalTk_GlbNotTkMu","m_mumu_GlobalTk_GlbNotTkMu",500,0,200);
  m_mumu_GlobalTk_GlbNotTkMu_ss = fs->make<TH1F>("m_mumu_GlobalTk_GlbNotTkMu_ss","m_mumu_GlobalTk_GlbNotTkMu_ss",500,0,200);

  m_mutrk_GlobalTk = fs->make<TH1F>("m_mutrk_GlobalTk","m_mutrk_GlobalTk",500,0,200);
  m_mutrk_GlobalTk_ss = fs->make<TH1F>("m_mutrk_GlobalTk_ss","m_mutrk_GlobalTk_ss",500,0,200);

  m_mumu_os_effcor = fs->make<TH1F>("m_mumu_os_effcor","m_mumu_os_effcor",500,0,200);
  m_mumu_ss_effcor = fs->make<TH1F>("m_mumu_ss_effcor","m_mumu_ss_effcor",500,0,200);
  m_mumu_2Global_ss_effcor = fs->make<TH1F>("m_mumu_2Global_ss_effcor","m_mumu_2Global_ss_effcor",500,0,200);
  m_mumu_2Global_effcor = fs->make<TH1F>("m_mumu_2Global_effcor","m_mumu_2Global_effcor",500,0,200);
  m_mutrk_GlobalTk_ss_effcor = fs->make<TH1F>("m_mutrk_GlobalTk_ss_effcor","m_mutrk_GlobalTk_ss_effcor",500,0,200);
  m_mutrk_GlobalTk_effcor = fs->make<TH1F>("m_mutrk_GlobalTk_effcor","m_mutrk_GlobalTk_effcor",500,0,200);

  eta_mu = fs->make<TH1F>("eta_mu","eta_mu",50,-2.5,2.5);
  eta_mu_mutrk = fs->make<TH1F>("eta_mu_mutrk","eta_mu_mutrk",50,-2.5,2.5);
  pt_mu_mutrk = fs->make<TH1F>("pt_mu_mutrk","pt_mu_mutrk",100,0,100);
  m_mumu_effcor = fs->make<TH1F>("m_mumu_effcor","m_mumu_effcor",500,0,200);
  pt_mu_effcor = fs->make<TH1F>("pt_mu_effcor","pt_mu_effcor",100,0,100);
  eta_mu_effcor = fs->make<TH1F>("eta_mu_effcor","eta_mu_effcor",50,-2.5,2.5);
  pt_low_mu = fs->make<TH1F>("pt_low_mu","pt_low_mu",20,0,20);
  pt_low_mu_effcor = fs->make<TH1F>("pt_low_mu_effcor","pt_low_mu_effcor",20,0,20);

  // mu pt
  pt_mu_Gen = fs->make<TH1F>("pt_mu_Gen","pt_mu_Gen",100,0,100);

  pt_mu = fs->make<TH1F>("pt_mu","pt_mu",100,0,100);
  pt_mu_Glb = fs->make<TH1F>("pt_mu_Glb","pt_mu_Glb",100,0,100);
  pt_mu_Glb_TkMu = fs->make<TH1F>("pt_mu_Glb_TkMu","pt_mu_Glb_TkMu",100,0,100);
  pt_mu_Glb_SA = fs->make<TH1F>("pt_mu_Glb_SA","pt_mu_Glb_SA",100,0,100);
  pt_mu_TkMu = fs->make<TH1F>("pt_mu_TkMu","pt_mu_TkMu",100,0,100);
  pt_mu_TkMu_Glb = fs->make<TH1F>("pt_mu_TkMu_Glb","pt_mu_TkMu_Glb",100,0,100);
  pt_mu_TkMu_SA = fs->make<TH1F>("pt_mu_TkMu_SA","pt_mu_TkMu_SA",100,0,100);
  pt_mu_SA = fs->make<TH1F>("pt_mu_SA","pt_mu_SA",100,0,100);
  pt_mu_SA_Glb = fs->make<TH1F>("pt_mu_SA_Glb","pt_mu_SA_Glb",100,0,100);
  pt_mu_SA_TkMu = fs->make<TH1F>("pt_mu_SA_TkMu","pt_mu_SA_TkMu",100,0,100);

  // PAT pt mu
  pt_patmu = fs->make<TH1F>("pt_patmu","pt_patmu",100,0,100);
  pt_patmu_Glb = fs->make<TH1F>("pt_patmu_Glb","pt_patmu_Glb",100,0,100);
  pt_patmu_Glb_TkMu = fs->make<TH1F>("pt_patmu_Glb_TkMu","pt_patmu_Glb_TkMu",100,0,100);
  pt_patmu_GlbNotTkMu = fs->make<TH1F>("pt_patmu_GlbNotTkMu","pt_patmu_GlbNotTkMu",100,0,100);
  pt_patmu_Glb_SA = fs->make<TH1F>("pt_patmu_Glb_SA","pt_patmu_Glb_SA",100,0,100);
  pt_patmu_TkMu = fs->make<TH1F>("pt_patmu_TkMu","pt_patmu_TkMu",100,0,100);
  pt_patmu_TkMu_Glb = fs->make<TH1F>("pt_patmu_TkMu_Glb","pt_patmu_TkMu_Glb",100,0,100);
  pt_patmu_TkMuNotGlb = fs->make<TH1F>("pt_patmu_TkMuNotGlb","pt_patmu_TkMuNotGlb",100,0,100);
  pt_patmu_TkMu_SA = fs->make<TH1F>("pt_patmu_TkMu_SA","pt_patmu_TkMu_SA",100,0,100);
  pt_patmu_SA = fs->make<TH1F>("pt_patmu_SA","pt_patmu_SA",100,0,100);
  pt_patmu_SA_Glb = fs->make<TH1F>("pt_patmu_SA_Glb","pt_patmu_SA_Glb",100,0,100);
  pt_patmu_SA_TkMu = fs->make<TH1F>("pt_patmu_SA_TkMu","pt_patmu_SA_TkMu",100,0,100);

  ptgen_patmuMC_Glb = fs->make<TH1F>("ptgen_patmuMC_Glb","ptgen_patmuMC_Glb",100,0,100);
  ptgen_patmuMC_Glb_TkMu = fs->make<TH1F>("ptgen_patmuMC_Glb_TkMu","ptgen_patmuMC_Glb_TkMu",100,0,100);
  ptgen_patmuMC_Glb_SA = fs->make<TH1F>("ptgen_patmuMC_Glb_SA","ptgen_patmuMC_Glb_SA",100,0,100);
  ptgen_patmuMC_GlbNotTkMu = fs->make<TH1F>("ptgen_patmuMC_GlbNotTkMu","ptgen_patmuMC_GlbNotTkMu",100,0,100);
  ptgen_patmuMC_TkMu = fs->make<TH1F>("ptgen_patmuMC_TkMu","ptgen_patmuMC_TkMu",100,0,100);
  ptgen_patmuMC_TkMu_SA = fs->make<TH1F>("ptgen_patmuMC_TkMu_SA","ptgen_patmuMC_TkMu_SA",100,0,100);
  ptgen_patmuMC_TkMuNotGlb = fs->make<TH1F>("ptgen_patmuMC_TkMuNotGlb","ptgen_patmuMC_TkMuNotGlb",100,0,100);
  ptgen_patmuMC_SA = fs->make<TH1F>("ptgen_patmuMC_SA","ptgen_patmuMC_SA",100,0,100);


  // pt reco - pt gen 
  Deltapt_muMC = fs->make<TH1F>("Deltapt_muMC","Deltapt_muMC",50,-5,5);
  Deltapt_patmuMC = fs->make<TH1F>("Deltapt_patmuMC","Deltapt_patmuMC",50,-5,5);
  Deltapt_patmuMC_Glb = fs->make<TH1F>("Deltapt_patmuMC_Glb","Deltapt_patmuMC_Glb",50,-5,5);
  Deltapt_patmuMC_Glb_TkMu = fs->make<TH1F>("Deltapt_patmuMC_Glb_TkMu","Deltapt_patmuMC_Glb_TkMu",50,-5,5);
  Deltapt_patmuMC_GlbNotTkMu = fs->make<TH1F>("Deltapt_patmuMC_GlbNotTkMu","Deltapt_patmuMC_GlbNotTkMu",50,-5,5);
  Deltapt_patmuMC_TkMu = fs->make<TH1F>("Deltapt_patmuMC_TkMu","Deltapt_patmuMC_TkMu",50,-5,5);
  Deltapt_patmuMC_TkMuNotGlb = fs->make<TH1F>("Deltapt_patmuMC_TkMuNotGlb","Deltapt_patmuMC_TkMuNotGlb",50,-5,5);
  Deltapt_patmuMC_SA = fs->make<TH1F>("Deltapt_patmuMC_SA","Deltapt_patmuMC_SA",100,-10,10);
  Deltapt_patmuMC_Glb_SA = fs->make<TH1F>("Deltapt_patmuMC_Glb_SA","Deltapt_patmuMC_Glb_SA",50,-5,5);
  Deltapt_patmuMC_TkMu_SA = fs->make<TH1F>("Deltapt_patmuMC_TkMu_SA","Deltapt_patmuMC_TkMu_SA",50,-5,5);

  // eta vs pt 
  etavspt_mu_Gen = fs->make<TH2F>("etavspt_mu_Gen","etavspt_mu_Gen",100,0,100,50,-2.5,2.5);
  etavspt_mu = fs->make<TH2F>("etavspt_mu","etavspt_mu",100,0,100,50,-2.5,2.5);
  etavspt_muNoMC = fs->make<TH2F>("etavspt_muNoMC","etavspt_muNoMC",100,0,100,50,-2.5,2.5);
  etavspt_mu_Glb = fs->make<TH2F>("etavspt_mu_Glb","etavspt_mu_Glb",100,0,100,50,-2.5,2.5);
  etavspt_mu_TkMu = fs->make<TH2F>("etavspt_mu_TkMu","etavspt_mu_TkMu",100,0,100,50,-2.5,2.5);
  etavspt_mu_SA = fs->make<TH2F>("etavspt_mu_SA","etavspt_mu_SA",100,0,100,50,-2.5,2.5);
  etavspt_patmu = fs->make<TH2F>("etavspt_patmu","etavspt_patmu",100,0,100,50,-2.5,2.5);
  etavspt_patmu_Glb = fs->make<TH2F>("etavspt_patmu_Glb","etavspt_patmu_Glb",100,0,100,50,-2.5,2.5); 
  etavspt_patmu_TkMu = fs->make<TH2F>("etavspt_patmu_TkMu","etavspt_patmu_TkMu",100,0,100,50,-2.5,2.5);
  etavspt_patmu_SA = fs->make<TH2F>("etavspt_patmu_SA","etavspt_patmu_SA",100,0,100,50,-2.5,2.5);
etavspt_patmuNoMC = fs->make<TH2F>("etavspt_patmuNoMC","etavspt_patmuNoMC",100,0,100,50,-2.5,2.5);
etavspt_patmuNoMC_Glb = fs->make<TH2F>("etavspt_patmuNoMC_Glb","etavspt_patmuNoMC_Glb",100,0,100,50,-2.5,2.5);
etavspt_patmuNoMC_TkMu = fs->make<TH2F>("etavspt_patmuNoMC_TkMu","etavspt_patmuNoMC_TkMu",100,0,100,50,-2.5,2.5);
etavspt_patmuNoMC_SA = fs->make<TH2F>("etavspt_patmuNoMC_SA","etavspt_patmuNoMC_SA",100,0,100,50,-2.5,2.5);

  etavspt_GlobalTkMu_lowm = fs->make<TH2F>("etavspt_GlobalTkMu_lowm","etavspt_GlobalTkMu_lowm",100,0,100,50,-2.5,2.5);
  etavspt_GlobalTkMu_lowm_ss = fs->make<TH2F>("etavspt_GlobalTkMu_lowm_ss","etavspt_GlobalTkMu_lowm_ss",100,0,100,50,-2.5,2.5);
  etavspt_GlobalTkMuSub = fs->make<TH2F>("etavspt_GlobalTkMuSub","etavspt_GlobalTkMuSub",100,0,100,50,-2.5,2.5);
  etavspt_GlobalTkMuSub_ss = fs->make<TH2F>("etavspt_GlobalTkMuSub_ss","etavspt_GlobalTkMuSub_ss",100,0,100,50,-2.5,2.5);
  etavspt_GlobalTk_GlbNotTkMu = fs->make<TH2F>("etavspt_GlobalTk_GlbNotTkMu","etavspt_GlobalTk_GlbNotTkMu",100,0,100,50,-2.5,2.5);
  etavspt_GlobalTk_GlbNotTkMu_ss = fs->make<TH2F>("etavspt_GlobalTk_GlbNotTkMu_ss","etavspt_GlobalTk_GlbNotTkMu_ss",100,0,100,50,-2.5,2.5); 
  etavspt_GlobalTk_TkMuNotGlb = fs->make<TH2F>("etavspt_GlobalTk_TkMuNotGlb","etavspt_GlobalTk_TkMuNotGlb",100,0,100,50,-2.5,2.5);
  etavspt_GlobalTk_TkMuNotGlb_ss = fs->make<TH2F>("etavspt_GlobalTk_TkMuNotGlb_ss","etavspt_GlobalTk_TkMuNotGlb_ss",100,0,100,50,-2.5,2.5);

  // charge*pt 
  qpt_mu_Gen = fs->make<TH1F>("qpt_mu_Gen","qpt_mu_Gen",200,-100,100);
  qpt_patmu = fs->make<TH1F>("qpt_patmu","qpt_patmu",200,-100,100);
  qpt_patmu_Glb = fs->make<TH1F>("qpt_patmu_Glb","qpt_patmu_Glb",200,-100,100);
  qpt_patmu_TkMu = fs->make<TH1F>("qpt_patmu_TkMu","qpt_patmu_TkMu",200,-100,100);
  qpt_patmu_SA = fs->make<TH1F>("qpt_patmu_SA","qpt_patmu_SA",200,-100,100);
  qpt_mu = fs->make<TH1F>("qpt_mu","qpt_mu",200,-100,100);
  qpt_mu_Glb = fs->make<TH1F>("qpt_mu_Glb","qpt_mu_Glb",200,-100,100);
  qpt_mu_TkMu = fs->make<TH1F>("qpt_mu_TkMu","qpt_mu_TkMu",200,-100,100);
  qpt_mu_SA = fs->make<TH1F>("qpt_mu_SA","qpt_mu_SA",200,-100,100);
  
 
//
  m_mumu_2mni = fs->make<TH1F>("m_mumu_2mni","m_mumu_2mni",500,0,200);
  m_mutrk_2mni = fs->make<TH1F>("m_mutrk_2mni","m_mutrk_2mni",500,0,200);
  m_musamu_2mni = fs->make<TH1F>("m_musamu_2mni","m_musamu_2mni",500,0,200);
                                                                                                                                                        
  eta_mu_2mni = fs->make<TH1F>("eta_mu_2mni","eta_mu_2mni",50,-2.5,2.5);
  pt_mu_2mni = fs->make<TH1F>("pt_mu_2mni","pt_mu_2mni",100,0,100);
  eta_mu_mutrk_2mni = fs->make<TH1F>("eta_mu_mutrk_2mni","eta_mu_mutrk_2mni",50,-2.5,2.5);
  pt_mu_mutrk_2mni = fs->make<TH1F>("pt_mu_mutrk_2mni","pt_mu_mutrk_2mni",100,0,100);
  m_mumu_effcor_2mni = fs->make<TH1F>("m_mumu_effcor_2mni","m_mumu_effcor_2mni",500,0,200);
  pt_mu_effcor_2mni = fs->make<TH1F>("pt_mu_effcor_2mni","pt_mu_effcor_2mni",100,0,100);
  eta_mu_effcor_2mni = fs->make<TH1F>("eta_mu_effcor_2mni","eta_mu_effcor_2mni",50,-2.5,2.5);
  pt_low_mu_2mni = fs->make<TH1F>("pt_low_mu_2mni","pt_low_mu_2mni",20,0,20);
  pt_low_mu_effcor_2mni = fs->make<TH1F>("pt_low_mu_effcor_2mni","pt_low_mu_effcor_2mni",20,0,20);


                                                                                                                            
  triggerbits = fs->make<TH1F>("triggerbits","triggerbits",25,0,25);

  triggerbits->GetXaxis()->SetBinLabel(1,"HLT_L1Mu"); 
  triggerbits->GetXaxis()->SetBinLabel(2,"HLT_L1MuOpen");
  triggerbits->GetXaxis()->SetBinLabel(3,"HLT_L2Mu9");
  triggerbits->GetXaxis()->SetBinLabel(4,"HLT_IsoMu9");
  triggerbits->GetXaxis()->SetBinLabel(5,"HLT_IsoMu11");
  triggerbits->GetXaxis()->SetBinLabel(6,"HLT_IsoMu13");
  triggerbits->GetXaxis()->SetBinLabel(7,"HLT_IsoMu15");
  triggerbits->GetXaxis()->SetBinLabel(8,"HLT_Mu3");
  triggerbits->GetXaxis()->SetBinLabel(9,"HLT_Mu5");
  triggerbits->GetXaxis()->SetBinLabel(10,"HLT_Mu7");
 triggerbits->GetXaxis()->SetBinLabel(11,"HLT_Mu9");
  triggerbits->GetXaxis()->SetBinLabel(12,"HLT_Mu11");
  triggerbits->GetXaxis()->SetBinLabel(13,"HLT_Mu13");
  triggerbits->GetXaxis()->SetBinLabel(14,"HLT_Mu15");
  triggerbits->GetXaxis()->SetBinLabel(15,"HLT_Mu15_L1Mu7");
  triggerbits->GetXaxis()->SetBinLabel(16,"HLT_Mu15_Vtx2cm");
  triggerbits->GetXaxis()->SetBinLabel(17,"HLT_Mu15_Vtx2mm");
  triggerbits->GetXaxis()->SetBinLabel(18,"HLT_DoubleIsoMu3");
  triggerbits->GetXaxis()->SetBinLabel(19,"HLT_DoubleMu3");
  triggerbits->GetXaxis()->SetBinLabel(20,"HLT_DoubleMu3_BJPsi");
  triggerbits->GetXaxis()->SetBinLabel(21,"HLT_DoubleMu4_BJPsi");
  triggerbits->GetXaxis()->SetBinLabel(22,"HLT_TripleMu3_TauTo3Mu");
  triggerbits->GetXaxis()->SetBinLabel(23,"HLT_Mu14_Jet50");
  triggerbits->GetXaxis()->SetBinLabel(24,"HLT_Mu5_TripleJet30");


}


DimuonsAnalyzer::~DimuonsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if ( ! SAFileName_.empty() ) {
     theefffile->Close();
   }
   if ( ! TKFileName_.empty() ) {
     theSAefffile->Close();
   }
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
DimuonsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Service<TFileService> fs;


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif


  // this analyzer produces a small root file with basic candidates and some MC information
  // some additional print statements
  nEvt++;
  TwoMuNonIso = false ;

  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
    std::cout << "reading event " << nEvt << std::endl;
 
  // HLT path
  edm::Handle<edm::TriggerResults> hltResults ;
  iEvent.getByLabel(InputTag("TriggerResults::HLT"),hltResults) ;
  trigNames.init(*hltResults) ;

  for (unsigned int i=0; i<trigNames.size(); i++)
    {
      //cout << "=== trigNames.triggerNames " << trigNames.triggerNames().at(i) <<endl;
      if ( trigNames.triggerNames().at(i) == "HLT_L1Mu" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(0);
        }
       if ( trigNames.triggerNames().at(i) == "HLT_L1MuOpen" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(1);
        }
       if ( trigNames.triggerNames().at(i) == "HLT_L2Mu9" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(2);
        }
       if ( trigNames.triggerNames().at(i) == "HLT_IsoMu9" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(3);
        }
       if ( trigNames.triggerNames().at(i) == "HLT_IsoMu11" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(4);
        }
       if ( trigNames.triggerNames().at(i) == "HLT_IsoMu13" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(5);
        }
       if ( trigNames.triggerNames().at(i) == "HLT_IsoMu15" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(6);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu3" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(7);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu5" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(8);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu7" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(9);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu9" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(10);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu11" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(11);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu13" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(12);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu15" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(13);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu15_L1Mu7" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(14);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu15_Vtx2cm" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(15);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu15_Vtx2mm" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(16);
        }

      if ( trigNames.triggerNames().at(i) == "HLT_DoubleIsoMu3" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(17);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_DoubleMu3" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(18);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_DoubleMu3_BJPsi" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(19);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_DoubleMu4_BJPsi" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(20);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_TripleMu3_TauTo3Mu" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(21);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu14_Jet50" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(22);
        }
      if ( trigNames.triggerNames().at(i) == "HLT_Mu5_TripleJet30" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(23);
        }

    }
       
   // generator info
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);
   int nGen = genParticles->size();
   int ngenCand=0;
   double DimuonGen_mumu_mass=0;
     
   for (int i=0; i<nGen; i++) {
    const Candidate &genCand = (*genParticles)[i];
    if ( genCand.status() == 3 && genCand.numberOfDaughters() == 3 ){
          ngenCand++;
          const Candidate * dauGen0 = genCand.daughter(0);
          const Candidate * dauGen1 = genCand.daughter(1);
          const Candidate * dauGen2 = genCand.daughter(2);
          Particle::LorentzVector pZ(0, 0, 0, 0);
	  Particle::LorentzVector muplusp4 = getParticleP4(-13,dauGen0,dauGen1,dauGen2);
	  Particle::LorentzVector muminusp4 = getParticleP4(13,dauGen0,dauGen1,dauGen2);
	  pZ = muplusp4 + muminusp4;
          DimuonGen_mumu_mass = pZ.mass();  // dimuons mass at generator level
          CandpdgId->Fill(genCand.pdgId());
          if (genCand.pdgId() == 23 ) {
            m_massGenZ->Fill(genCand.mass());
            m_mumu_GenZ->Fill(DimuonGen_mumu_mass);
          } 
         // check for di-muons
         if ( ( abs(dauGen0->pdgId()) == 13 && abs(dauGen0->pdgId()) == 13 )  || ( abs(dauGen0->pdgId()) == 13 && abs(dauGen2->pdgId()) == 13 )  || ( abs(dauGen1->pdgId()) == 13 && abs(dauGen2->pdgId()) == 13 )  )  {
            cout << " genCand ID "<< genCand.pdgId() << ", dauGen0=" << dauGen0->pdgId() << ", dauGen1=" << dauGen1->pdgId() << ", dauGen2=" << dauGen2->pdgId() << endl;

           m_massGen->Fill(genCand.mass());
           m_mumu_Gen->Fill(DimuonGen_mumu_mass);
           pt_mu_Gen->Fill(muplusp4.pt());
           pt_mu_Gen->Fill(muminusp4.pt());
           qpt_mu_Gen->Fill(muplusp4.pt());
           qpt_mu_Gen->Fill(muminusp4.pt()*(-1));
           etavspt_mu_Gen->Fill( muplusp4.pt(), -log(tan(muplusp4.theta()/2.)) );
           etavspt_mu_Gen->Fill( muminusp4.pt(), -log(tan(muminusp4.theta()/2.)) );

         }

    } 

   }


  // MC map
  Handle<GenParticleMatch> muMatchMap;
  iEvent.getByLabel("allDimuonsMCMatch",muMatchMap);
  const GenParticleMatch & mcMatch = * muMatchMap;

// PAT Layer1 Mu
  Handle<edm::View<pat::Muon> > PATmuons;
  iEvent.getByLabel("selectedLayer1Muons",PATmuons);
  edm::View<pat::Muon>::const_iterator PATmuon;
  nPATMuon=0;
  nPATMuonMCMatch=0;
  TrackRef TkTrack;
  TrackRef SATrack;
  size_t nPATmu = PATmuons->size();
  cout << "nPATmu " << nPATmu <<endl;

  for( size_t i = 0; i < nPATmu  && nPATMuon<MUONMAX; i++ ) {
     //const Candidate & PATmuon = (*PATmuons)[i];
     //ok cout << "i PatMu pt "<< PATmuon.pt() << " TrackMu? " << PATmuon.isTrackerMuon() <<endl;
      RefToBase<pat::Muon> PATmuon = PATmuons->refAt(i);
      CandidateBaseRef muCandRef(PATmuon);
      GenParticleRef muGen= mcMatch[muCandRef];
      if( muGen.isNonnull() ) {
            nPATMuonMCMatch++;
            //cout << "PAT MATCH found!!! " << endl;
            double ptGen = muGen->pt();
            double ptReco = PATmuon->pt();
            cout << "PAT MATCH found: ptGen "<< ptGen << " ptReco " << ptReco << endl;
            Deltapt_patmuMC->Fill(PATmuon->pt()-muGen->pt());
      } else {
           etavspt_patmuNoMC->Fill(PATmuon->pt(),PATmuon->eta());
      }
        pt_patmu->Fill(PATmuon->pt());
        qpt_patmu->Fill(PATmuon->pt()*PATmuon->charge());
        etavspt_patmu->Fill(PATmuon->pt(),PATmuon->eta());

// CaloMuons: there are no caloMuon in this collection
        if ( PATmuon->isCaloMuon()){ cout << "PAT CaloMuon pt:" << PATmuon->pt() << endl;}
// Global
        if ( PATmuon->isGlobalMuon()){
            pt_patmu_Glb->Fill(PATmuon->pt());
            qpt_patmu_Glb->Fill(PATmuon->pt()*PATmuon->charge());
            etavspt_patmu_Glb->Fill(PATmuon->pt(),PATmuon->eta());
            if( muGen.isNonnull() ) {
               ptgen_patmuMC_Glb->Fill(muGen->pt());
               Deltapt_patmuMC_Glb->Fill(PATmuon->pt()-muGen->pt());
            } else {
              etavspt_patmuNoMC_Glb->Fill(PATmuon->pt(),PATmuon->eta());
            }
            if ( PATmuon->isTrackerMuon())    { // global & TkMu
                 pt_patmu_Glb_TkMu->Fill(PATmuon->pt());
                 if( muGen.isNonnull() ) {
                    ptgen_patmuMC_Glb_TkMu->Fill(muGen->pt()); 
                    Deltapt_patmuMC_Glb_TkMu->Fill(PATmuon->pt()-muGen->pt());
                 }
            } else {  // global & ! TkMu
                 pt_patmu_GlbNotTkMu->Fill(PATmuon->pt());
                 if( muGen.isNonnull() ) {
                    ptgen_patmuMC_GlbNotTkMu->Fill(muGen->pt()); 
                    Deltapt_patmuMC_GlbNotTkMu->Fill(PATmuon->pt()-muGen->pt());
                 }
            }
            if ( PATmuon->isStandAloneMuon()) { // global & SA
              pt_patmu_Glb_SA->Fill(PATmuon->pt());
              if( muGen.isNonnull() ) {
                  ptgen_patmuMC_Glb_SA->Fill(muGen->pt());
                  Deltapt_patmuMC_Glb_SA->Fill(PATmuon->pt()-muGen->pt());
              }
            } 
        }
// TrackerMu
        if ( PATmuon->isTrackerMuon()){
         // if ( PATmuon->track()->isNonnull()){ // add validity check
         TkTrack = PATmuon->track();
            pt_patmu_TkMu->Fill(TkTrack->pt());
            qpt_patmu_TkMu->Fill(TkTrack->pt()*TkTrack->charge());
            etavspt_patmu_TkMu->Fill(TkTrack->pt(),TkTrack->eta());
            if( muGen.isNonnull() ) {
               ptgen_patmuMC_TkMu->Fill(muGen->pt());
               Deltapt_patmuMC_TkMu->Fill(PATmuon->pt()-muGen->pt());
            } else {
              etavspt_patmuNoMC_TkMu->Fill(PATmuon->pt(),PATmuon->eta());
            }
            if ( PATmuon->isGlobalMuon())     { // TkMu & global
               pt_patmu_TkMu_Glb->Fill(TkTrack->pt());
            } else { // TkMu & ! global
               pt_patmu_TkMuNotGlb->Fill(TkTrack->pt());
               if( muGen.isNonnull() ) {
                    ptgen_patmuMC_TkMuNotGlb->Fill(muGen->pt());
                    Deltapt_patmuMC_TkMuNotGlb->Fill(PATmuon->pt()-muGen->pt());
               }
            } 
            if ( PATmuon->isStandAloneMuon()) {
              pt_patmu_TkMu_SA->Fill(TkTrack->pt());
              if( muGen.isNonnull() ) {
                  ptgen_patmuMC_TkMu_SA->Fill(muGen->pt());
                  Deltapt_patmuMC_TkMu_SA->Fill(PATmuon->pt()-muGen->pt());
              }

             }
         //}
        }
// SA
        if ( PATmuon->isStandAloneMuon()){
            SATrack = PATmuon->standAloneMuon();
            pt_patmu_SA->Fill(SATrack->pt());
            qpt_patmu_SA->Fill(SATrack->pt()*SATrack->charge());
            etavspt_patmu_SA->Fill(SATrack->pt(),SATrack->eta());
            if( muGen.isNonnull() ) {
               ptgen_patmuMC_SA->Fill(muGen->pt());
               Deltapt_patmuMC_SA->Fill(SATrack->pt()-muGen->pt());
            } else {
              etavspt_patmuNoMC_SA->Fill(SATrack->pt(),SATrack->eta());
            }
            if ( PATmuon->isGlobalMuon())  pt_patmu_SA_Glb->Fill(SATrack->pt());
            if ( PATmuon->isTrackerMuon()) pt_patmu_SA_TkMu->Fill(SATrack->pt());
        }

   nPATMuon++;
  }

// Reco Mu
// reco::Muon* muon = ...
//if (muon->track()->isNonnull()) {
// reco::TrackRef innerTrackRef = muon->track()
// reco::TrackRef outerTrackRef = muon->standAloneTrack()

   // dimuons skim
   Handle<reco::CompositeCandidateCollection> dimuons;
   iEvent.getByLabel("dimuons",dimuons);
   //cout << "-- dimuons size : " << dimuons->size() << endl;

   reco::CompositeCandidateCollection::const_iterator dimuon;
   nDimuonCand=0;
   nMatchedMu=0;
   TrackRef SATrack1;
   TrackRef SATrack2;
//   TrackRef Track1;
//   TrackRef Track2;

   for( dimuon = dimuons->begin(); dimuon != dimuons->end() && nDimuonCand<DIMUONMAX; ++ dimuon ) {

     DimuonCand_goodmumu_mass[nDimuonCand]=dimuon->mass();
     MuonCand_pt1[nDimuonCand] = dimuon->daughter(0)->pt();
     MuonCand_pt2[nDimuonCand] = dimuon->daughter(1)->pt();
     MuonCand_eta1[nDimuonCand] = dimuon->daughter(0)->eta();
     MuonCand_eta2[nDimuonCand] = dimuon->daughter(1)->eta();


     bool SSAccept = false;

     bool GlobalsAccept = false;
     bool GlobalSAAccept = false;
     bool SAsAccept = false;
     bool GlobalTrackerMuonAccept = false;
     bool TrackerMuonsAccept = false;
     bool GlobalCaloMuonAccept = false; 
     bool GlobalTrackerMuonSubsetAccept = false;
     bool GlobalTrackerMuonGlbNotTkMu = false;
     bool GlobalTrackerMuonTkMuNotGlb = false;
     bool Cand1TkMuNotGlb = false;
     bool Cand2TkMuNotGlb = false;
     bool Cand2GlbNotTkMu = false;
     bool Cand1GlbNotTkMu = false;

     double effTk1 = 1.0;
     double effTk2 = 1.0;
     double effSA1 = 1.0; 
     double effSA2 = 1.0;
     // read SA TnP efficiency 
     if ( ! SAFileName_.empty() ) {
         effSA1 = readefficiency(dimuon->daughter(0)->pt(),dimuon->daughter(0)->eta());
         effSA2 = readefficiency(dimuon->daughter(1)->pt(),dimuon->daughter(1)->eta());
     }
     // read TnP efficiency
     if ( ! TKFileName_.empty() ) {
         effTk1 = readTKefficiency(dimuon->daughter(0)->pt(),dimuon->daughter(0)->eta());
         effTk2 = readTKefficiency(dimuon->daughter(1)->pt(),dimuon->daughter(1)->eta());
     }

     double mu1weight = effSA1 * effTk1 ;
     double mu2weight = effSA2 * effTk2 ;
     double dimuweight = (1.0/(mu1weight*mu2weight));


  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100)) {
     cout << "candidate with mass " << dimuon->mass() <<endl;  
     cout << "Weight = " << (1.0/mu1weight) << ", " << (1.0/mu2weight) << ", " << dimuweight << endl;
     cout << "isGlobalMuon = " << dimuon->daughter(0)->isGlobalMuon() << ", " << dimuon->daughter(1)->isGlobalMuon()<< endl;     
     cout << "isTrackerMuon = " << dimuon->daughter(0)->isTrackerMuon() << ", " << dimuon->daughter(1)->isTrackerMuon()<< endl;
   }
     // selection based on mu types    
// same sign
     if ( dimuon->daughter(0)->charge() == dimuon->daughter(1)->charge() ) SSAccept = true;
// global + global
     if ( dimuon->daughter(0)->isGlobalMuon() && dimuon->daughter(1)->isGlobalMuon() ) GlobalsAccept = true;
// global + trackerMu
     if ( ( dimuon->daughter(0)->isTrackerMuon() && dimuon->daughter(1)->isGlobalMuon() ) || ( dimuon->daughter(1)->isTrackerMuon() && dimuon->daughter(0)->isGlobalMuon() )) GlobalTrackerMuonAccept = true;
// global && !trackerMu + trackerMu && ! global
     if ( dimuon->daughter(0)->isTrackerMuon() && !dimuon->daughter(0)->isGlobalMuon() ) Cand1TkMuNotGlb=true;
     if ( dimuon->daughter(1)->isTrackerMuon() && !dimuon->daughter(1)->isGlobalMuon() ) Cand2TkMuNotGlb=true;
     if ( dimuon->daughter(1)->isGlobalMuon() && !dimuon->daughter(1)->isTrackerMuon() ) Cand2GlbNotTkMu=true;
     if ( dimuon->daughter(0)->isGlobalMuon() && !dimuon->daughter(0)->isTrackerMuon() ) Cand1GlbNotTkMu=true;
     if ( ( Cand1TkMuNotGlb && Cand2GlbNotTkMu ) || (Cand2TkMuNotGlb && Cand1GlbNotTkMu) ) GlobalTrackerMuonSubsetAccept=true;
     if ( ( Cand1GlbNotTkMu && dimuon->daughter(1)->isTrackerMuon() ) || ( Cand2GlbNotTkMu && dimuon->daughter(0)->isTrackerMuon() ) ) GlobalTrackerMuonGlbNotTkMu=true;
     if ( ( Cand1TkMuNotGlb && dimuon->daughter(1)->isGlobalMuon() ) || ( Cand2TkMuNotGlb && dimuon->daughter(0)->isGlobalMuon() ) ) GlobalTrackerMuonTkMuNotGlb=true;
// trackerMu + trackerMu
     if ( dimuon->daughter(0)->isTrackerMuon() && dimuon->daughter(1)->isTrackerMuon() ) TrackerMuonsAccept = true;
// global + caloMu
     if ( ( dimuon->daughter(0)->isCaloMuon() && dimuon->daughter(1)->isGlobalMuon() ) || ( dimuon->daughter(1)->isCaloMuon() && dimuon->daughter(0)->isGlobalMuon() )) GlobalCaloMuonAccept = true;
// global + SA
     if ( ( dimuon->daughter(0)->isStandAloneMuon() && dimuon->daughter(1)->isGlobalMuon() ) || ( dimuon->daughter(1)->isStandAloneMuon() && dimuon->daughter(0)->isGlobalMuon() )) GlobalSAAccept = true;
// SA + SA
     if ( dimuon->daughter(0)->isStandAloneMuon() && dimuon->daughter(1)->isStandAloneMuon() ) SAsAccept = true;

     int charge1=0;
     int charge2=0;

//
// MC match
//
     const reco::Candidate * dau0 = dimuon->daughter(0);
     reco::CandidateBaseRef mu1 = dau0->masterClone();
     GenParticleRef muGen1= mcMatch[mu1];
     if( muGen1.isNonnull() ) {
            nMatchedMu++;
            //cout << "MATCH found!!! " << endl;
            //cout << "pt gen "<< muGen1->pt() << "ptReco " <<  MuonCand_pt1[nDimuonCand] <<endl;
            Deltapt_muMC->Fill(MuonCand_pt1[nDimuonCand]-muGen1->pt());
     } else {
       etavspt_muNoMC->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
     }
                                                                                                                       
     const reco::Candidate * dau1 = dimuon->daughter(1);
     reco::CandidateBaseRef mu2 = dau1->masterClone();
     GenParticleRef muGen2= mcMatch[mu2];
     if( muGen2.isNonnull() ) {
            nMatchedMu++;
            //cout << "MATCH found!!! " << endl;
            //cout << "pt gen "<< muGen2->pt() << "ptReco" << MuonCand_pt2[nDimuonCand] <<endl;
            Deltapt_muMC->Fill(MuonCand_pt2[nDimuonCand]-muGen2->pt()); 
     } else {
       etavspt_muNoMC->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
     }


// global + SA : with SA momentum (that double the statistic for global+global events)
     double Dimuon_SAG_mass=0.;
     if ( GlobalSAAccept ) { 
      Dimuon_SAG_mass = 0.;
      if ( dimuon->daughter(0)->isStandAloneMuon() )
      {
          SATrack1 = dimuon->daughter(0)->get<TrackRef,reco::StandAloneMuonTag>();
         //SATrack1 = dimuon->daughter(0)->standAloneMuon();
         if (! SATrack1.isNull()){ 
          double cand_p = dimuon->daughter(1)->p() + SATrack1->p();
          double cand_px = dimuon->daughter(1)->px() + SATrack1->px();
          double cand_py = dimuon->daughter(1)->py() + SATrack1->py();
          double cand_pz = dimuon->daughter(1)->pz() + SATrack1->pz();
          Dimuon_SAG_mass = cand_p*cand_p - cand_px*cand_px - cand_py*cand_py - cand_pz*cand_pz;
          Dimuon_SAG_mass = (Dimuon_SAG_mass>0) ? sqrt(Dimuon_SAG_mass) : 0.;
          charge2 = dimuon->daughter(1)->charge();
          charge1 = SATrack1->charge();
          //cout << " Dimuon_SAG_mass " <<Dimuon_SAG_mass <<endl;
          if ( charge1 == charge2 ) {
            m_mumu_GlobalSA_ss->Fill(Dimuon_SAG_mass);
          } else {
            m_mumu_GlobalSA->Fill(Dimuon_SAG_mass);
          }
         }
      } 
      if ( dimuon->daughter(1)->isStandAloneMuon() ) {
        SATrack2 = dimuon->daughter(1)->get<TrackRef,reco::StandAloneMuonTag>(); 
        if (! SATrack2.isNull()){ 
          double cand_p = dimuon->daughter(0)->p() + SATrack2->p();
          double cand_px = dimuon->daughter(0)->px() + SATrack2->px();
          double cand_py = dimuon->daughter(0)->py() + SATrack2->py();
          double cand_pz = dimuon->daughter(0)->pz() + SATrack2->pz();
          Dimuon_SAG_mass = cand_p*cand_p - cand_px*cand_px - cand_py*cand_py - cand_pz*cand_pz;
          Dimuon_SAG_mass = (Dimuon_SAG_mass>0) ? sqrt(Dimuon_SAG_mass) : 0.;
          charge1 = dimuon->daughter(0)->charge();
          charge2 = SATrack2->charge();
          //cout << " Dimuon_SAG_mass " <<Dimuon_SAG_mass <<endl;
          if ( charge1 == charge2 ) {
            m_mumu_GlobalSA_ss->Fill(Dimuon_SAG_mass);
          } else {
            m_mumu_GlobalSA->Fill(Dimuon_SAG_mass);
          }
        } 
      }
     }

// SA + SA
    double Dimuon_2SA_mass=0.;
    if ( SAsAccept ) {
         SATrack1 = dimuon->daughter(0)->get<TrackRef,reco::StandAloneMuonTag>();
         SATrack2 = dimuon->daughter(1)->get<TrackRef,reco::StandAloneMuonTag>();
         if ( (! SATrack1.isNull()) && (! SATrack2.isNull()) ){
          double cand_p = SATrack1->p() + SATrack2->p();
          double cand_px = SATrack1->px() + SATrack2->px();
          double cand_py = SATrack1->py() + SATrack2->py();
          double cand_pz = SATrack1->pz() + SATrack2->pz();
          Dimuon_2SA_mass = cand_p*cand_p - cand_px*cand_px - cand_py*cand_py - cand_pz*cand_pz;
          Dimuon_2SA_mass = (Dimuon_2SA_mass>0) ? sqrt(Dimuon_2SA_mass) : 0.;
          if ( SATrack1->charge() == SATrack2->charge() ) {
            m_mumu_2SA_ss->Fill(Dimuon_2SA_mass); 
          } else { 
            m_mumu_2SA->Fill(Dimuon_2SA_mass);
          }
         }
    }

// global + trackerMuon with trackerMuon momentum from track()
     double Dimuon_trkG_mass=0.;
     if ( GlobalTrackerMuonAccept ) {
       
       if ( dimuon->daughter(0)->isTrackerMuon() ) {
         // Track1 = dimuon->daughter(0)->track();
         TrackRef Track1 = dimuon->daughter(0)->get<TrackRef>();
         if ( Track1.isNonnull() ) {
          cout << "daughter(0) TrackerMuon:  daughter pt="<<dimuon->daughter(0)->pt() << " track pt="<<Track1->pt()<<endl;
          double cand_p = dimuon->daughter(1)->p() + Track1->p();
          double cand_px = dimuon->daughter(1)->px() + Track1->px();
          double cand_py = dimuon->daughter(1)->py() + Track1->py();
          double cand_pz = dimuon->daughter(1)->pz() + Track1->pz();
          Dimuon_trkG_mass = cand_p*cand_p - cand_px*cand_px - cand_py*cand_py - cand_pz*cand_pz;
          Dimuon_trkG_mass = (Dimuon_trkG_mass>0) ? sqrt(Dimuon_trkG_mass) : 0.;
          charge2 = dimuon->daughter(1)->charge();
          charge1 = Track1->charge();
          //cout << " Dimuon_trkG_mass " <<Dimuon_trkG_mass <<endl;
           if ( charge1 == charge2 ) {
             m_mumu_GlobalTkMu_trk_ss->Fill(Dimuon_trkG_mass);
           } else {
             m_mumu_GlobalTkMu_trk->Fill(Dimuon_trkG_mass);
           }
         }
       }
       if ( dimuon->daughter(1)->isTrackerMuon() ) {
         TrackRef Track2 = dimuon->daughter(1)->get<TrackRef>();
         if ( Track2.isNonnull() ) {
          cout << "daughter(1) TrackerMuon:  daughter pt="<<dimuon->daughter(1)->pt() << " track pt="<<Track2->pt()<<endl;
          double cand_p = dimuon->daughter(0)->p() + Track2->p();
          double cand_px = dimuon->daughter(0)->px() + Track2->px();
          double cand_py = dimuon->daughter(0)->py() + Track2->py();
          double cand_pz = dimuon->daughter(0)->pz() + Track2->pz();
          Dimuon_trkG_mass = cand_p*cand_p - cand_px*cand_px - cand_py*cand_py - cand_pz*cand_pz;
          Dimuon_trkG_mass = (Dimuon_trkG_mass>0) ? sqrt(Dimuon_trkG_mass) : 0.;
          charge2 = dimuon->daughter(0)->charge();
          charge1 = Track2->charge();
          //cout << " Dimuon_trkG_mass " <<Dimuon_trkG_mass <<endl;
           if ( charge1 == charge2 ) {
             m_mumu_GlobalTkMu_trk_ss->Fill(Dimuon_trkG_mass);
           } else {
             m_mumu_GlobalTkMu_trk->Fill(Dimuon_trkG_mass);
           }
         }
       }
     }
// trackerMuon + trackerMuon with trackerMuon momentum from track()
     double Dimuon_2trk_mass=0.;
     if ( TrackerMuonsAccept ){
         TrackRef tTrack1 = dimuon->daughter(0)->get<TrackRef>();
         TrackRef tTrack2 = dimuon->daughter(1)->get<TrackRef>();

        if ( tTrack1.isNonnull() && tTrack2.isNonnull() ) {
          double cand_p = tTrack2->p() + tTrack1->p();
          double cand_px = tTrack2->px() + tTrack1->px();
          double cand_py = tTrack2->py() + tTrack1->py();
          double cand_pz = tTrack2->pz() + tTrack1->pz();
          Dimuon_2trk_mass = cand_p*cand_p - cand_px*cand_px - cand_py*cand_py - cand_pz*cand_pz;
          Dimuon_2trk_mass = (Dimuon_2trk_mass>0) ? sqrt(Dimuon_2trk_mass) : 0.;
          charge2 = tTrack2->charge();
          charge1 = tTrack1->charge();
          //cout << " Dimuon_2trk_mass " <<Dimuon_2trk_mass <<endl;
          if ( charge1 == charge2 ) {
            m_mumu_2TkMu_trk_ss->Fill(Dimuon_2trk_mass);
          } else {
            m_mumu_2TkMu_trk->Fill(Dimuon_2trk_mass);
          }
        }
     }


     if(mu1weight == 0 || mu2weight == 0)
       {
         mu1weight = 1.0;
         mu2weight = 1.0;
         dimuweight = 1.0;
       }
                                                                                                                               
     pt_mu_effcor->Fill(MuonCand_pt1[nDimuonCand],1.0/mu1weight);
     pt_mu_effcor->Fill(MuonCand_pt2[nDimuonCand],1.0/mu2weight);
     eta_mu_effcor->Fill(MuonCand_eta1[nDimuonCand],1.0/mu1weight);
     eta_mu_effcor->Fill(MuonCand_eta2[nDimuonCand],1.0/mu2weight);
                                                                                                                               
     if ( ! SSAccept ) {
       m_mumu_os->Fill(DimuonCand_goodmumu_mass[nDimuonCand]); 
       m_mumu_os_effcor->Fill(DimuonCand_goodmumu_mass[nDimuonCand],dimuweight);
       if ( GlobalsAccept ) {
        m_mumu_2Global->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
        m_mumu_2Global_effcor->Fill(DimuonCand_goodmumu_mass[nDimuonCand],dimuweight);
       }
       // if ( GlobalSAAccept ) m_mumu_GlobalSA->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
       //if ( GlobalSAAccept ) m_mumu_GlobalSA->Fill(Dimuon_SAG_mass);
       if ( GlobalTrackerMuonAccept ) {
          m_mumu_GlobalTkMu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]); 
          if ( TrackerMuonsAccept ) m_mumu_GlobalTkMuAND2TkMu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
          if ( GlobalsAccept ) m_mumu_GlobalTkMuAND2Glb->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
          if ( GlobalTrackerMuonSubsetAccept ){
           m_mumu_GlobalTkSub->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
           etavspt_GlobalTkMuSub->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTkMuSub->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
          }
          if ( GlobalTrackerMuonGlbNotTkMu ) {
           m_mumu_GlobalTk_GlbNotTkMu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);            
           etavspt_GlobalTk_GlbNotTkMu->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTk_GlbNotTkMu->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
          }
          if ( GlobalTrackerMuonTkMuNotGlb ) {
           m_mumu_GlobalTk_TkMuNotGlb->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
           etavspt_GlobalTk_TkMuNotGlb->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTk_TkMuNotGlb->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
          }

          if (DimuonCand_goodmumu_mass[nDimuonCand]<2.5){
           etavspt_GlobalTkMu_lowm->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTkMu_lowm->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
          }
       }
       if ( TrackerMuonsAccept )  m_mumu_2TkMu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
       if ( GlobalCaloMuonAccept ) m_mumu_GlobalCaloMu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
     } else {
       m_mumu_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
       m_mumu_ss_effcor->Fill(DimuonCand_goodmumu_mass[nDimuonCand],dimuweight);
       if ( GlobalsAccept ) {
        m_mumu_2Global_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
        m_mumu_2Global_ss_effcor->Fill(DimuonCand_goodmumu_mass[nDimuonCand],dimuweight);
       }
       //if ( GlobalSAAccept ) m_mumu_GlobalSA_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
       //if ( GlobalSAAccept ) m_mumu_GlobalSA_ss->Fill(Dimuon_SAG_mass);
       if ( GlobalTrackerMuonAccept ){
         m_mumu_GlobalTkMu_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);        
         if ( TrackerMuonsAccept ) m_mumu_GlobalTkMuAND2TkMu_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
         if ( GlobalsAccept ) m_mumu_GlobalTkMuAND2Glb_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
         if ( GlobalTrackerMuonSubsetAccept ) {
           m_mumu_GlobalTkSub_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
           etavspt_GlobalTkMuSub_ss->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTkMuSub_ss->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
         }

          if ( GlobalTrackerMuonGlbNotTkMu ) {
           m_mumu_GlobalTk_GlbNotTkMu_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
           etavspt_GlobalTk_GlbNotTkMu_ss->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTk_GlbNotTkMu_ss->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
          }
          if ( GlobalTrackerMuonTkMuNotGlb ) {
           m_mumu_GlobalTk_TkMuNotGlb_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
           etavspt_GlobalTk_TkMuNotGlb_ss->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTk_TkMuNotGlb_ss->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
          }

         if (DimuonCand_goodmumu_mass[nDimuonCand]<2.5){
           etavspt_GlobalTkMu_lowm_ss->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
           etavspt_GlobalTkMu_lowm_ss->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
          }

       }
       if ( TrackerMuonsAccept )  m_mumu_2TkMu_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
       if ( GlobalCaloMuonAccept ) m_mumu_GlobalCaloMu_ss->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
     }


     m_mumu_effcor->Fill(DimuonCand_goodmumu_mass[nDimuonCand],dimuweight); 
     m_mumu->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);

// muon pt

     pt_mu->Fill(MuonCand_pt1[nDimuonCand]);
     pt_mu->Fill(MuonCand_pt2[nDimuonCand]);
     etavspt_mu->Fill(MuonCand_pt1[nDimuonCand],MuonCand_eta1[nDimuonCand]);
     etavspt_mu->Fill(MuonCand_pt2[nDimuonCand],MuonCand_eta2[nDimuonCand]);
     qpt_mu->Fill(MuonCand_pt1[nDimuonCand]*dimuon->daughter(0)->charge());
     qpt_mu->Fill(MuonCand_pt2[nDimuonCand]*dimuon->daughter(1)->charge());

//MuonCand1
     if ( dimuon->daughter(0)->isGlobalMuon()){
            pt_mu_Glb->Fill(dimuon->daughter(0)->pt());
            qpt_mu_Glb->Fill(dimuon->daughter(0)->pt()*dimuon->daughter(0)->charge());
            etavspt_mu_Glb->Fill(dimuon->daughter(0)->pt(),dimuon->daughter(0)->eta());
            if ( dimuon->daughter(0)->isTrackerMuon())    pt_mu_Glb_TkMu->Fill(dimuon->daughter(0)->pt());
            if ( dimuon->daughter(0)->isStandAloneMuon()) pt_mu_Glb_SA->Fill(dimuon->daughter(0)->pt());
     }
     if ( dimuon->daughter(0)->isTrackerMuon()){
            TrackRef TrackCand1 = dimuon->daughter(0)->get<TrackRef>();
         if ( TrackCand1.isNonnull() ) {
            pt_mu_TkMu->Fill(TrackCand1->pt());
            qpt_mu_TkMu->Fill(TrackCand1->pt()*TrackCand1->charge());
            etavspt_mu_TkMu->Fill(TrackCand1->pt(),TrackCand1->eta());
            if ( dimuon->daughter(0)->isGlobalMuon())     pt_mu_TkMu_Glb->Fill(TrackCand1->pt());
            if ( dimuon->daughter(0)->isStandAloneMuon()) pt_mu_TkMu_SA->Fill(TrackCand1->pt());
         }
     }
     if ( dimuon->daughter(0)->isStandAloneMuon()){
            TrackRef SACand1 = dimuon->daughter(0)->get<TrackRef,reco::StandAloneMuonTag>();
         if (! SACand1.isNull()){
            pt_mu_SA->Fill(SACand1->pt());
            qpt_mu_SA->Fill(SACand1->pt()*SACand1->charge());
            etavspt_mu_SA->Fill(SACand1->pt(),SACand1->eta());
            if ( dimuon->daughter(0)->isGlobalMuon())  pt_mu_SA_Glb->Fill(SACand1->pt());
            if ( dimuon->daughter(0)->isTrackerMuon()) pt_mu_SA_TkMu->Fill(SACand1->pt());
         }
     }
//MuonCand2
     if ( dimuon->daughter(1)->isGlobalMuon()){
            pt_mu_Glb->Fill(dimuon->daughter(1)->pt());
            qpt_mu_Glb->Fill(dimuon->daughter(1)->pt()*dimuon->daughter(1)->charge());
            etavspt_mu_Glb->Fill(dimuon->daughter(1)->pt(),dimuon->daughter(1)->eta());
            if ( dimuon->daughter(1)->isTrackerMuon())    pt_mu_Glb_TkMu->Fill(dimuon->daughter(1)->pt());
            if ( dimuon->daughter(1)->isStandAloneMuon()) pt_mu_Glb_SA->Fill(dimuon->daughter(1)->pt());
     }
     if ( dimuon->daughter(1)->isTrackerMuon()){
            TrackRef TrackCand2 = dimuon->daughter(1)->get<TrackRef>();
          if ( TrackCand2.isNonnull() ) { 
            pt_mu_TkMu->Fill(TrackCand2->pt());
            qpt_mu_TkMu->Fill(TrackCand2->pt()*TrackCand2->charge());
            etavspt_mu_TkMu->Fill(TrackCand2->pt(),TrackCand2->eta());
            if ( dimuon->daughter(1)->isGlobalMuon())     pt_mu_TkMu_Glb->Fill(TrackCand2->pt());
            if ( dimuon->daughter(1)->isStandAloneMuon()) pt_mu_TkMu_SA->Fill(TrackCand2->pt());
          }
     }
     if ( dimuon->daughter(1)->isStandAloneMuon()){
            TrackRef SACand2 = dimuon->daughter(1)->get<TrackRef,reco::StandAloneMuonTag>();
         if (! SACand2.isNull()) {
            pt_mu_SA->Fill(SACand2->pt());
            qpt_mu_SA->Fill(SACand2->pt()*SACand2->charge());
            etavspt_mu_SA->Fill(SACand2->pt(),SACand2->eta());
            if ( dimuon->daughter(1)->isGlobalMuon())  pt_mu_SA_Glb->Fill(SACand2->pt());
            if ( dimuon->daughter(1)->isTrackerMuon()) pt_mu_SA_TkMu->Fill(SACand2->pt());
         }
     }


     eta_mu->Fill(MuonCand_eta1[nDimuonCand]);
     eta_mu->Fill(MuonCand_eta2[nDimuonCand]);
                                                                                                                               
     pt_low_mu->Fill(MuonCand_pt1[nDimuonCand]);
     pt_low_mu->Fill(MuonCand_pt2[nDimuonCand]);
     pt_low_mu_effcor->Fill(MuonCand_pt1[nDimuonCand],1.0/mu1weight);
     pt_low_mu_effcor->Fill(MuonCand_pt2[nDimuonCand],1.0/mu2weight);
          
     if ( TwoMuNonIso )
     {
     pt_mu_effcor_2mni->Fill(MuonCand_pt1[nDimuonCand],1.0/mu1weight);
     pt_mu_effcor_2mni->Fill(MuonCand_pt2[nDimuonCand],1.0/mu2weight);
     eta_mu_effcor_2mni->Fill(MuonCand_eta1[nDimuonCand],1.0/mu1weight);
     eta_mu_effcor_2mni->Fill(MuonCand_eta2[nDimuonCand],1.0/mu2weight);                                                                                                                                                   
     m_mumu_effcor_2mni->Fill(DimuonCand_goodmumu_mass[nDimuonCand],dimuweight);
     m_mumu_2mni->Fill(DimuonCand_goodmumu_mass[nDimuonCand]);
     pt_mu_2mni->Fill(MuonCand_pt1[nDimuonCand]);
     pt_mu_2mni->Fill(MuonCand_pt2[nDimuonCand]);
     eta_mu_2mni->Fill(MuonCand_eta1[nDimuonCand]);
     eta_mu_2mni->Fill(MuonCand_eta2[nDimuonCand]);
                                                                                                                                                        
     pt_low_mu_2mni->Fill(MuonCand_pt1[nDimuonCand]);
     pt_low_mu_2mni->Fill(MuonCand_pt2[nDimuonCand]);
     pt_low_mu_effcor_2mni->Fill(MuonCand_pt1[nDimuonCand],1.0/mu1weight);
     pt_low_mu_effcor_2mni->Fill(MuonCand_pt2[nDimuonCand],1.0/mu2weight);

     } 
                                                                                                                     
     nDimuonCand++;
   }


   // dimuonsOneTrack skim : 
   Handle<reco::CompositeCandidateCollection> trkdimuons;
   iEvent.getByLabel("dimuonsOneTrack",trkdimuons);
   //cout << "-- dimuonsOneTrack size : " << trkdimuons->size() << endl;
   reco::CompositeCandidateCollection::const_iterator trkdimuon;
   nTrkDimuonCand=0;
   for( trkdimuon = trkdimuons->begin(); trkdimuon != trkdimuons->end() && nTrkDimuonCand<DIMUONMAX; ++ trkdimuon ) {

     bool trkSSAccept = false; 
     bool GlobalTrackerAccept = false;
     DimuonCand_mumuonetrack_mass[nTrkDimuonCand]=trkdimuon->mass();

     if ( trkdimuon->daughter(0)->charge() == trkdimuon->daughter(1)->charge() ) trkSSAccept = true;
     if ( trkdimuon->daughter(0)->isGlobalMuon() || trkdimuon->daughter(1)->isGlobalMuon() ) GlobalTrackerAccept = true;    

     m_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);

     double effSAa = 1.0;
     double effSAb = 1.0;
     double effTka = 1.0;
     double effTkb = 1.0;

     // read SA TnP efficiency
     if ( ! SAFileName_.empty() ) {
       if ( trkdimuon->daughter(0)->isGlobalMuon() ){
        effSAa = readefficiency(trkdimuon->daughter(0)->pt(),trkdimuon->daughter(0)->eta());
       }
       if ( trkdimuon->daughter(1)->isGlobalMuon() ){
        effSAb = readefficiency(trkdimuon->daughter(1)->pt(),trkdimuon->daughter(1)->eta());
       }
     }
     // read TK TnP efficiency
     if ( ! TKFileName_.empty() ) {
         effTka = readTKefficiency(trkdimuon->daughter(0)->pt(),trkdimuon->daughter(0)->eta());
         effTkb = readTKefficiency(trkdimuon->daughter(1)->pt(),trkdimuon->daughter(1)->eta());
     }

     double muaweight = effSAa * effTka ;
     double mubweight = effSAb * effTkb ;
     double dimutkweight = (1.0/(muaweight*mubweight));


     if ( ! trkSSAccept ) {
        if ( GlobalTrackerAccept ) {
          m_mutrk_GlobalTk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
          m_mutrk_GlobalTk_effcor->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand],dimutkweight);
        }
     } else {
        if ( GlobalTrackerAccept ) {
         m_mutrk_GlobalTk_ss->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
         m_mutrk_GlobalTk_ss_effcor->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand],dimutkweight);
        }
     }
          if ( TwoMuNonIso )
     {
     m_mutrk_2mni->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
     } 
                                                                                                                         
     nTrkDimuonCand++;
   }


}


Particle::LorentzVector DimuonsAnalyzer::getParticleP4(const int ipart, const Candidate * dauGen0, const Candidate * dauGen1, const Candidate * dauGen2) 
{
  int partId0 = dauGen0->pdgId();
  int partId1 = dauGen1->pdgId();
  int partId2 = dauGen2->pdgId();
  Particle::LorentzVector p4part(0.,0.,0.,0.);
  if (partId0 == ipart) {
    for(size_t k = 0; k < dauGen0->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen0->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	p4part = dauMuGen->p4();
      }
    }
  }
  if (partId1 == ipart) {
    for(size_t k = 0; k < dauGen1->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen1->daughter(k);
      if(dauMuGen->pdgId() == ipart && dauMuGen->status() ==1) {
	p4part = dauMuGen->p4();
      }
    }
  }
  if (partId2 == ipart) {
    for(size_t k = 0; k < dauGen2->numberOfDaughters(); ++k) {
      const Candidate * dauMuGen = dauGen2->daughter(k);
      if(abs(dauMuGen->pdgId()) == ipart && dauMuGen->status() ==1) {
	p4part = dauMuGen->p4();
      }
    }
  }
  return p4part;
}

/*
bool DimuonsAnalyzer::isZmumu(const Candidate *genCand )
{
 bool Zmumu=false;
 if((genCand.pdgId() == 23)&&(genCand.status() == 3)) {
 const Candidate * dauGen0 = genCand.daughter(0);
 const Candidate * dauGen1 = genCand.daughter(1);
 const Candidate * dauGen2 = genCand.daughter(2);
                                                                                                                      
  int partId0 = dauGen0->pdgId();
  int partId1 = dauGen1->pdgId();
  int partId2 = dauGen2->pdgId();
  bool muplusFound=false;
  bool muminusFound=false;
  bool ZFound=false;
  if (partId0==13 || partId1==13 || partId2==13) muminusFound=true;
  if (partId0==-13 || partId1==-13 || partId2==-13) muplusFound=true;
  if (partId0==23 || partId1==23 || partId2==23) ZFound=true;
  Zmumu=muplusFound*muminusFound*ZFound;
 }
 return Zmumu                                                                                                                    
}*/

double DimuonsAnalyzer::readefficiency2D(double pt, double eta)
{
  double theeff = 1.0;
  theSAefffile->cd();

  TH2F *hptetaeff = new TH2F();
  hptetaeff = (TH2F *)theSAefffile->Get("sbs_eff_Pt_Eta");


  for(int i = 0;i < hptetaeff->GetNbinsX(); i++)
    {
      for(int j = 0;j < hptetaeff->GetNbinsY(); j++)
        {
          if(pt > hptetaeff->GetXaxis()->GetBinLowEdge(i+1) &&
             pt < (hptetaeff->GetXaxis()->GetBinLowEdge(i+1) + hptetaeff->GetXaxis()->GetBinWidth(i+1)))
            {
              if(eta > hptetaeff->GetYaxis()->GetBinLowEdge(j+1) &&
                 eta < (hptetaeff->GetYaxis()->GetBinLowEdge(j+1) + hptetaeff->GetYaxis()->GetBinWidth(j+1)))
                {
                  theeff = hptetaeff->GetBinContent(i+1,j+1);
                }
            }
        }
    }

  return theeff;
}


double DimuonsAnalyzer::readefficiency(double pt, double eta)
{
  double theeff = 1.0;
  theSAefffile->cd();

  TH1F *hpteff = new TH1F();
  if ( ( sampletype_ == "Upsilon") || ( sampletype_ == "JPsi") || ( sampletype_ == "MuonPT5")){
    hpteff = (TH1F *)theSAefffile->Get("sbs_eff_Pt");

    for(int i = 0;i < hpteff->GetNbinsX(); i++)
    {
      if(pt > hpteff->GetXaxis()->GetBinLowEdge(i+1) &&
             pt < (hpteff->GetXaxis()->GetBinLowEdge(i+1) + hpteff->GetXaxis()->GetBinWidth(i+1)))
            {
              theeff = hpteff->GetBinContent(i+1);
            }
    } 
  } else {
    hpteff = (TH1F *)theSAefffile->Get("fit_eff_Pt");
   //double binhigh = 0.0;
    for(int i = 0;i < hpteff->GetNbinsX(); i++)
    {
      if(  pt > hpteff->GetXaxis()->GetBinLowEdge(i+1) &&
           pt < (hpteff->GetXaxis()->GetBinLowEdge(i+1) + hpteff->GetXaxis()->GetBinWidth(i+1)) )
            {
              theeff = hpteff->GetBinContent(i+1);
            }
    }
  } 
  return theeff;
 
}

double DimuonsAnalyzer::readTKefficiency(double pt, double eta)
{
  double theeff = 1.0;
                                                                                
  //This is testing the readback of Tag&Probe efficiencies
  theefffile->cd();
  TH1F *hpteff = new TH1F();
  // use sideband subtracted because fit is screwed for pt >90Gev
  //hpteff = (TH1F *)theefffile->Get("fit_eff_Pt");
  hpteff = (TH1F *)theefffile->Get("sbs_eff_Pt");
  for(int i = 0;i < hpteff->GetNbinsX(); i++)
  {
      if(  pt > hpteff->GetXaxis()->GetBinLowEdge(i+1) &&
           pt < (hpteff->GetXaxis()->GetBinLowEdge(i+1) + hpteff->GetXaxis()->GetBinWidth(i+1)) )
            {
              theeff = hpteff->GetBinContent(i+1);
            }
  }
  return theeff;
}


// ------------ method called once each job just before starting event loop  ------------
void 
DimuonsAnalyzer::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DimuonsAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DimuonsAnalyzer);
