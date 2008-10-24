//; -*- C++ -*-
//
//; Package:    DimuonsZselAnalyzer
// Class:      DimuonsZselAnalyzer
// 
/**\class DimuonsZselAnalyzer DimuonsZselAnalyzer.cc Dimuons/DimuonsZselAnalyzer/src/DimuonsZselAnalyzer.cc

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
 
#include "DataFormats/Common/interface/ValueMap.h"

// TFileService files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
                                                                                                                              
#include <TFile.h>
#include <TH1D.h>
#include <TObject.h>
#include <TTree.h>
                                                                                                                               
using namespace std;
using namespace reco;
using namespace edm;
using namespace HepMC;

//
// class decleration
//

class DimuonsZselAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DimuonsZselAnalyzer(const edm::ParameterSet&);
      ~DimuonsZselAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      double readefficiency(double pt, double eta);
      double readTKefficiency(double pt, double eta);

      // ----------member data ---------------------------
  int nEvt;// used to count the number of events
                                                                                                                               
  TFile *theefffile;
  TFile *theSAefffile;                                                                                                                             
  // to be used for root output tree
  TH1F *m_mumu;
  TH1F *m_mutrk;
  TH1F *m_musamu;

  TH1F *m_mumu_effcor;

  TH1F *pt_mu;
  TH1F *pt_mu_effcor;
  TH1F *pt_low_mu;
  TH1F *pt_low_mu_effcor;
  TH1F *eta_mu;
  TH1F *eta_mu_effcor;
  TH1F *pt_mu_mutrk;
  TH1F *eta_mu_mutrk;
  
  TH1F *pt_mumu;
  TH1F *pt_mumu_effcor; 
  TH1F *eta_mumu;
  TH1F *eta_mumu_effcor;
  TH1F *y_mumu;
  TH1F *y_mumu_effcor;

  TH1F *triggerbits;

  TH1F *zsel_m_mumu;
  TH1F *zsel_m_mutrk;
  TH1F *zsel_m_musamu;                                                                                                       
  TH1F *zsel_m_mumu_effcor;
  TH1F *zsel_m_mumu_geffcor;                                                                                                         
  TH1F *zsel_pt_mu;
  TH1F *zsel_pt_mu_effcor;
  TH1F *zsel_pt_low_mu;
  TH1F *zsel_pt_low_mu_effcor;
  TH1F *zsel_eta_mu;
  TH1F *zsel_eta_mu_effcor;
  TH1F *zsel_pt_mu_mutrk;
  TH1F *zsel_eta_mu_mutrk;
                                                                                                         
  TH1F *zsel_pt_mumu;
  TH1F *zsel_pt_mumu_effcor;
  TH1F *zsel_pt_mumu_geffcor;
  TH1F *zsel_eta_mumu;
  TH1F *zsel_eta_mumu_effcor;
  TH1F *zsel_eta_mumu_geffcor;
  TH1F *zsel_y_mumu;
  TH1F *zsel_y_mumu_effcor;
  TH1F *zsel_y_mumu_geffcor;

                                                                                                                               
  edm::TriggerNames trigNames ;
 
  int nDimuonCand, nTrkDimuonCand, nSADimuonCand;
  int nTrig;
/*
  int DIMUONMAX;// used to set maximum of arrays
  double DimuonCand_goodmumu_mass[10];
  double DimuonCand_mumuonetrack_mass[10];
  double DimuonCand_mumuonesamuon_mass[10];
  double MuonCand_pt1[10], MuonCand_pt2[10];
  double MuonCand_eta1[10], MuonCand_eta2[10];

  double DimuonCand_goodmumu_pt[10];                                           
  double DimuonCand_goodmumu_eta[10];                               
  double DimuonCand_goodmumu_y[10];          
  */
  double DimuonCand_goodmumu_mass;   
  double DimuonCand_mumuonetrack_mass;
  double DimuonCand_mumuonesamuon_mass;
  double DimuonCand_goodmumu_pt;
  double DimuonCand_goodmumu_eta;
  double DimuonCand_goodmumu_y;

                                         
  int nbinspt;
  double ptlow;
  double pthigh;
  double binwidth;

  edm::InputTag muIso1_, muIso2_;
  double isocut_, etacut_, ptcut_,  minZmass_, maxZmass_;

  double SAEff ; 
  double TkEff ; 
  double IsoEff ; 
  double MuEff ; 
  
  string TKFileName_  ;
  string SAFileName_  ;

};

//
// constants, enums and typedefs
//

// typedef edm::AssociationVector<reco::CandidateRefProd, std::vector<double> > IsolationCollection;

  typedef edm::ValueMap<float> IsolationCollection; 
//

// static data member definitions
//

//
// constructors and destructor
//
DimuonsZselAnalyzer::DimuonsZselAnalyzer(const edm::ParameterSet& iConfig) :  
  muIso1_(iConfig.getParameter<edm::InputTag>("muonIsolations1")),
  muIso2_(iConfig.getParameter<edm::InputTag>("muonIsolations2")),
  isocut_(iConfig.getParameter<double>( "isocut" ) ),
  etacut_(iConfig.getParameter<double>( "etacut" ) ),
  ptcut_(iConfig.getParameter<double>( "ptcut" ) ),
  minZmass_(iConfig.getParameter<double>( "minZmass" )),
  maxZmass_(iConfig.getParameter<double>( "maxZmass" ))  {


  TKFileName_ = iConfig.getUntrackedParameter<string> ("TkEffFile");
  SAFileName_ = iConfig.getUntrackedParameter<string> ("SAEffFile");

  //now do what ever initialization is needed
  nEvt=0;
  // DIMUONMAX=10;
  nTrig = 4;
  nbinspt = 20; ptlow = 0.0; pthigh = 100.0; binwidth = (pthigh-ptlow)/nbinspt;

  SAEff =  0.937 ;
  TkEff = 0.9966 ;
  IsoEff = 0.974 ;
  MuEff = SAEff * TkEff ;

 
  theefffile = TFile::Open(TKFileName_.c_str());
  theSAefffile = TFile::Open(SAFileName_.c_str());

//  theefffile = TFile::Open("/bohome/fanfani/CRAB/trackingmuon_eff_testBin8.root");
//  theSAefffile = TFile::Open("/bohome/fanfani/CRAB/SAmuon_eff_testBin8.root");

  //theefffile = TFile::Open("muon_eff_test.root");
           
  //TFileService
  edm::Service<TFileService> fs;
                                                                                                                               
  m_mumu = fs->make<TH1F>("m_mumu","m_mumu",500,0,200);
  m_mutrk = fs->make<TH1F>("m_mutrk","m_mutrk",500,0,200);
  m_musamu = fs->make<TH1F>("m_musamu","m_musamu",500,0,200);
  // Z                   
  pt_mumu = fs->make<TH1F>("pt_mumu","pt_mumu",100,0,100);
  eta_mumu = fs->make<TH1F>("eta_mumu","eta_mumu",50,-2.5,2.5);
  y_mumu = fs->make<TH1F>("y_mumu","y_mumu",50,-5,5);
  pt_mumu_effcor = fs->make<TH1F>("pt_mumu_effcor","pt_mumu_effcor",100,0,100);
  eta_mumu_effcor = fs->make<TH1F>("eta_mumu_effcor","eta_mumu_effcor",50,-2.5,2.5);
  y_mumu_effcor = fs->make<TH1F>("y_mumu_effcor","y_mumu_effcor",50,-5,5);                                                                                              
  eta_mu = fs->make<TH1F>("eta_mu","eta_mu",50,-2.5,2.5);
  pt_mu = fs->make<TH1F>("pt_mu","pt_mu",100,0,100);
  eta_mu_mutrk = fs->make<TH1F>("eta_mu_mutrk","eta_mu_mutrk",50,-2.5,2.5);
  pt_mu_mutrk = fs->make<TH1F>("pt_mu_mutrk","pt_mu_mutrk",100,0,100);
  m_mumu_effcor = fs->make<TH1F>("m_mumu_effcor","m_mumu_effcor",500,0,200);
  pt_mu_effcor = fs->make<TH1F>("pt_mu_effcor","pt_mu_effcor",100,0,100);
  eta_mu_effcor = fs->make<TH1F>("eta_mu_effcor","eta_mu_effcor",50,-2.5,2.5);
  pt_low_mu = fs->make<TH1F>("pt_low_mu","pt_low_mu",20,0,20);
  pt_low_mu_effcor = fs->make<TH1F>("pt_low_mu_effcor","pt_low_mu_effcor",20,0,20);
   
/*
  triggerbits = fs->make<TH1F>("triggerbits","triggerbits",6,0,6);
  triggerbits->GetXaxis()->SetBinLabel(1,"HLT1MuonPrescalePt3");
  triggerbits->GetXaxis()->SetBinLabel(2,"HLT1MuonPrescalePt5");
  triggerbits->GetXaxis()->SetBinLabel(3,"HLT1MuonPrescalePt7x7");
  triggerbits->GetXaxis()->SetBinLabel(4,"HLT1MuonIso");
  triggerbits->GetXaxis()->SetBinLabel(5,"HLT1MuonNonIso15");
  triggerbits->GetXaxis()->SetBinLabel(6,"HLT2MuonNonIso");
*/

  zsel_m_mumu = fs->make<TH1F>("zsel_m_mumu","zsel_m_mumu",500,0,200);
  zsel_m_mutrk = fs->make<TH1F>("zsel_m_mutrk","zsel_m_mutrk",500,0,200);
  zsel_m_musamu = fs->make<TH1F>("zsel_m_musamu","zsel_m_musamu",500,0,200);
  // Z
  zsel_pt_mumu = fs->make<TH1F>("zsel_pt_mumu","zsel_pt_mumu",100,0,100);
  zsel_eta_mumu = fs->make<TH1F>("zsel_eta_mumu","zsel_eta_mumu",50,-2.5,2.5);
  zsel_y_mumu = fs->make<TH1F>("zsel_y_mumu","zsel_y_mumu",50,-5,5); 
  zsel_pt_mumu_effcor = fs->make<TH1F>("zsel_pt_mumu_effcor","zsel_pt_mumu_effcor",100,0,100);
  zsel_eta_mumu_effcor = fs->make < TH1F > ("zsel_eta_mumu_effcor", "zsel_eta_mumu_effcor", 50, -2.5, 2.5);
  zsel_y_mumu_effcor = fs->make < TH1F > ("zsel_y_mumu_effcor","zsel_y_mumu_effcor", 50, -5, 5);
  zsel_pt_mumu_geffcor = fs->make<TH1F>("zsel_pt_mumu_geffcor","zsel_pt_mumu_geffcor",100,0,100);
  zsel_eta_mumu_geffcor = fs->make < TH1F > ("zsel_eta_mumu_geffcor", "zsel_eta_mumu_geffcor", 50, -2.5, 2.5);
  zsel_y_mumu_geffcor = fs->make < TH1F > ("zsel_y_mumu_geffcor","zsel_y_mumu_geffcor", 50, -5, 5);

  zsel_eta_mu = fs->make<TH1F>("zsel_eta_mu","zsel_eta_mu",50,-2.5,2.5);
  zsel_pt_mu = fs->make<TH1F>("zsel_pt_mu","zsel_pt_mu",100,0,100);
  zsel_eta_mu_mutrk = fs->make<TH1F>("zsel_eta_mu_mutrk","zsel_eta_mu_mutrk",50,-2.5,2.5);
  zsel_pt_mu_mutrk = fs->make<TH1F>("zsel_pt_mu_mutrk","zsel_pt_mu_mutrk",100,0,100);
  zsel_m_mumu_effcor = fs->make<TH1F>("zsel_m_mumu_effcor","zsel_m_mumu_effcor",500,0,200);
  zsel_pt_mu_effcor = fs->make<TH1F>("zsel_pt_mu_effcor","zsel_pt_mu_effcor",100,0,100);
  zsel_eta_mu_effcor = fs->make<TH1F>("zsel_eta_mu_effcor","zsel_eta_mu_effcor",50,-2.5,2.5);
  zsel_m_mumu_geffcor = fs->make<TH1F>("zsel_m_mumu_geffcor","zsel_m_mumu_geffcor",500,0,200);
  zsel_pt_low_mu = fs->make<TH1F>("zsel_pt_low_mu","zsel_pt_low_mu",20,0,20);
  zsel_pt_low_mu_effcor = fs->make<TH1F>("zsel_pt_low_mu_effcor","zsel_pt_low_mu_effcor",20,0,20);

}


DimuonsZselAnalyzer::~DimuonsZselAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   theefffile->Close();
   theSAefffile->Close();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
DimuonsZselAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  if((nEvt%10==0 && nEvt<=100)||(nEvt%100==0 && nEvt>100))
    std::cout << "reading event " << nEvt << std::endl;
   
  //   std::cout << "Got Event" << std::endl;
 

/*
  edm::Handle<edm::TriggerResults> hltResults ;
  //  edm::Handle<l1extra::L1ParticleMapCollection> l1Results ;
 
  iEvent.getByLabel(InputTag("TriggerResults::HLT"),hltResults) ;
  //  iEvent.getByLabel("l1extraParticleMap",l1Results) ;
  trigNames.init(*hltResults) ;
 
  //  cout << "# of triggers = " << trigNames.size() << endl;
 
  for (unsigned int i=0; i<trigNames.size(); i++)
    {
      //      cout << "\tTrigger = " << trigNames.triggerNames().at(i) << endl;
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt3" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(0);
        }
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt5" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(1);
        }
      if ( trigNames.triggerNames().at(i) == "HLT1MuonPrescalePt7x7" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(2);
        }
      if ( trigNames.triggerNames().at(i) == "HLT1MuonIso" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(3);
        }
      if ( trigNames.triggerNames().at(i) == "HLT1MuonNonIso15" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(4);
        }
      if ( trigNames.triggerNames().at(i) == "HLT2MuonNonIso" )
        {
          if ( hltResults->accept(i) )
            triggerbits->Fill(5);
        }
    }
*/                                                                                                                               

/*
   Handle<reco::CompositeCandidateCollection> dimuons;
   iEvent.getByLabel("dimuons",dimuons);
*/
   Handle<edm::View<reco::Candidate> > zCands;
   iEvent.getByLabel("zToMuGlobal", zCands);
   size_t nZMuMu = zCands->size();

   reco::CompositeCandidateCollection::const_iterator dimuon;
// isolation
   Handle<IsolationCollection> hMuIso1_, hMuIso2_;
   iEvent.getByLabel(muIso1_, hMuIso1_);
   iEvent.getByLabel(muIso2_, hMuIso2_);

   nDimuonCand=0;
   int ZMuMu_selection = 0;
   bool ptcutAccept = false;
   bool etacutAccept = false;
   bool masscutAccept = false;
   bool isocutAccept = false;
   bool OSAccept = false;

/*
   for( dimuon = dimuons->begin(); dimuon != dimuons->end() && nDimuonCand<DIMUONMAX; ++ dimuon ) {
*/ 
   for( size_t i = 0; i < nZMuMu ; i++ ) {

   const Candidate & dimuon = (*zCands)[ i ];
 
     const reco::Candidate * dau0 = dimuon.daughter(0);
     const reco::Candidate * dau1 = dimuon.daughter(1);
/*
     const reco::Candidate * dau0 = dimuon->daughter(0);
     const reco::Candidate * dau1 = dimuon->daughter(1);
*/
     reco::CandidateBaseRef mu1 = dau0->masterClone();
     reco::CandidateBaseRef mu2 = dau1->masterClone();
     double iso1 = (*hMuIso1_)[mu1];
     double iso2 = (*hMuIso2_)[mu2];
/*
     double pt1m = dimuon->daughter(0)->pt();
     double pt2m = dimuon->daughter(1)->pt();
     double eta1m = dimuon->daughter(0)->eta();
     double eta2m = dimuon->daughter(1)->eta();
*/
     double pt1m = dimuon.daughter(0)->pt();
     double pt2m = dimuon.daughter(1)->pt();
     double eta1m = dimuon.daughter(0)->eta();
     double eta2m = dimuon.daughter(1)->eta();

     if ( pt1m > ptcut_ && pt2m > ptcut_ ) ptcutAccept = true;
     if ( fabs(eta1m)<etacut_ && fabs(eta2m)<etacut_ ) etacutAccept = true;
     if (iso1 < isocut_ && iso2 <isocut_) isocutAccept = true;
     if ( dimuon.mass() > minZmass_ && dimuon.mass() < maxZmass_ ) masscutAccept = true;
     if ( dimuon.daughter(0)->charge() != dimuon.daughter(1)->charge() ) OSAccept = true;   

     DimuonCand_goodmumu_mass=dimuon.mass();
     DimuonCand_goodmumu_pt=dimuon.pt();
     DimuonCand_goodmumu_eta=dimuon.eta();
     DimuonCand_goodmumu_y=dimuon.rapidity();
/*
     MuonCand_pt1[nDimuonCand] = pt1m;
     MuonCand_pt2[nDimuonCand] = pt2m;
     MuonCand_eta1[nDimuonCand] = eta1m;
     MuonCand_eta2[nDimuonCand] = eta2m;
*/                                                                                                                               

/*
     cout << "Candidate with mass = " << dimuon->mass() << endl;
     double mu1weight = readefficiency(dimuon->daughter(0)->pt(),dimuon->daughter(0)->eta());
     double mu2weight = readefficiency(dimuon->daughter(1)->pt(),dimuon->daughter(1)->eta());
*/
     cout << "Candidate with mass = " << dimuon.mass() << endl;

/*
     double effTk1 = readTKefficiency(dimuon.daughter(0)->pt(),dimuon.daughter(0)->eta());
     double effTk2 = readTKefficiency(dimuon.daughter(1)->pt(),dimuon.daughter(1)->eta());
     double effSA1 = readefficiency(dimuon.daughter(0)->pt(),dimuon.daughter(0)->eta());
     double effSA2 = readefficiency(dimuon.daughter(1)->pt(),dimuon.daughter(1)->eta());
*/
     double effTk1 = readTKefficiency(pt1m,eta1m);
     double effTk2 = readTKefficiency(pt2m,eta2m);
     double effSA1 = readefficiency(pt1m,eta1m);
     double effSA2 = readefficiency(pt2m,eta2m);
     double mu1weight = effSA1 * effTk1 ;
     double mu2weight = effSA2 * effTk2 ;

     double dimuweight = (1.0/(mu1weight*mu2weight));
     cout << "Weight = " << (1.0/mu1weight) << ", " << (1.0/mu2weight) << ", " << dimuweight << endl;
     cout << "ptcutAccept " << ptcutAccept << ", etacutAccept " << etacutAccept << ", masscutAccept " << masscutAccept << " isocutAccept " << isocutAccept << "OSAccept "<< OSAccept << endl; 
     cout << "pt1m=" << pt1m << " eta1m=" << eta1m << " effSA1=" << effSA1 << " effTk1=" << effTk1<< endl;               cout << "pt2m=" << pt2m << " eta2m=" << eta2m << " effSA2=" << effSA2 << " effTk2=" << effTk2<< endl;
     if(mu1weight == 0 || mu2weight == 0)
       {
         mu1weight = 1.0;
         mu2weight = 1.0;
         dimuweight = 1.0;
       }

     double gdimuweight = (1.0/(MuEff*MuEff)); 

     pt_mu_effcor->Fill(pt1m,1.0/mu1weight);
     pt_mu_effcor->Fill(pt2m,1.0/mu2weight);
     eta_mu_effcor->Fill(eta1m,1.0/mu1weight);
     eta_mu_effcor->Fill(eta2m,1.0/mu2weight);

     m_mumu_effcor->Fill(DimuonCand_goodmumu_mass,dimuweight);
     m_mumu->Fill(DimuonCand_goodmumu_mass);
     pt_mu->Fill(pt1m);
     pt_mu->Fill(pt2m);
     eta_mu->Fill(eta1m);
     eta_mu->Fill(eta2m);
     // Z candidate                                           
     pt_mumu_effcor->Fill(DimuonCand_goodmumu_pt,dimuweight);
     pt_mumu->Fill(DimuonCand_goodmumu_pt);
     eta_mumu_effcor->Fill(DimuonCand_goodmumu_eta,dimuweight);
     eta_mumu->Fill(DimuonCand_goodmumu_eta);                  
     y_mumu_effcor->Fill(DimuonCand_goodmumu_y,dimuweight);
     y_mumu->Fill(DimuonCand_goodmumu_y);
     //                                                           
     pt_low_mu->Fill(pt1m);
     pt_low_mu->Fill(pt2m);
     pt_low_mu_effcor->Fill(pt1m,1.0/mu1weight);
     pt_low_mu_effcor->Fill(pt2m,1.0/mu2weight);
 
    if ( ptcutAccept && etacutAccept && masscutAccept && isocutAccept && OSAccept ) {

     ZMuMu_selection++;
     zsel_pt_mu_effcor->Fill(pt1m,1.0/mu1weight);
     zsel_pt_mu_effcor->Fill(pt2m,1.0/mu2weight);
     zsel_eta_mu_effcor->Fill(eta1m,1.0/mu1weight);
     zsel_eta_mu_effcor->Fill(eta2m,1.0/mu2weight);
                                                                                                         
     zsel_m_mumu_effcor->Fill(DimuonCand_goodmumu_mass,dimuweight);
     zsel_m_mumu_geffcor->Fill(DimuonCand_goodmumu_mass,gdimuweight);
     zsel_m_mumu->Fill(DimuonCand_goodmumu_mass);
     zsel_pt_mu->Fill(pt1m);
     zsel_pt_mu->Fill(pt2m);
     zsel_eta_mu->Fill(eta1m);
     zsel_eta_mu->Fill(eta2m);
     // Z candidate
     zsel_pt_mumu_effcor->Fill(DimuonCand_goodmumu_pt,dimuweight);
     zsel_pt_mumu->Fill(DimuonCand_goodmumu_pt);
     zsel_eta_mumu_effcor->Fill(DimuonCand_goodmumu_eta,dimuweight);
     zsel_eta_mumu->Fill(DimuonCand_goodmumu_eta);
     zsel_y_mumu_effcor->Fill(DimuonCand_goodmumu_y,dimuweight);
     zsel_y_mumu->Fill(DimuonCand_goodmumu_y);

     zsel_pt_mumu_geffcor->Fill(DimuonCand_goodmumu_pt,gdimuweight);
     zsel_eta_mumu_geffcor->Fill(DimuonCand_goodmumu_eta,gdimuweight);     
     zsel_y_mumu_geffcor->Fill(DimuonCand_goodmumu_y,gdimuweight);
     //
     zsel_pt_low_mu->Fill(pt1m);
     zsel_pt_low_mu->Fill(pt2m);
     zsel_pt_low_mu_effcor->Fill(pt1m,1.0/mu1weight);
     zsel_pt_low_mu_effcor->Fill(pt2m,1.0/mu2weight);
     }


    nDimuonCand++;


   }

/*
   Handle<reco::CompositeCandidateCollection> trkdimuons;
   iEvent.getByLabel("dimuonsOneTrack",trkdimuons);
   reco::CompositeCandidateCollection::const_iterator trkdimuon;
   nTrkDimuonCand=0;
   for( trkdimuon = trkdimuons->begin(); trkdimuon != trkdimuons->end() ; ++ trkdimuon ) {
     DimuonCand_mumuonetrack_mass[nTrkDimuonCand]=trkdimuon->mass();
 
     m_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
     m_jpsi_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
     m_ups_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
     m_z_mutrk->Fill(DimuonCand_mumuonetrack_mass[nTrkDimuonCand]);
                                                                                                                               
     //     pt_mu_mutrk->Fill(dimuon->daughter(0)->pt());
     //     pt_mu_mutrk->Fill(dimuon->daughter(1)->pt());
     //     eta_mu_mutrk->Fill(dimuon->daughter(0)->eta());
     //     eta_mu_mutrk->Fill(dimuon->daughter(1)->eta());
                                                                                                                               
     nTrkDimuonCand++;
   }
*/


}

double DimuonsZselAnalyzer::readTKefficiency(double pt, double eta)
{
  double theeff = 1.0;
 
  //This is testing the readback of Tag&Probe efficiencies
  theefffile->cd();
  TH1F *hpteff = new TH1F();
  // use sideband subtracted because fit is screwed for pt >90Gev
  //hpteff = (TH1F *)theefffile->Get("fit_eff_Pt");
  hpteff = (TH1F *)theefffile->Get("sbs_eff_Pt");
  //  double binhigh = 0.0;
  double binlow = 10.0;
  double binwidth = 10. ;
  for(int i = 0;i < hpteff->GetNbinsX(); i++)
      {
        if(pt > (binlow + i * binwidth) && pt < ( binlow + (i+1)*binwidth))
      theeff = hpteff->GetBinContent(i+1);
      }
  return theeff;
}

double DimuonsZselAnalyzer::readefficiency(double pt, double eta)
{
  double theeff = 1.0;
  theSAefffile->cd();
  TH1F *hpteff = new TH1F();
  hpteff = (TH1F *)theSAefffile->Get("fit_eff_Pt");
  //double binhigh = 0.0;
  double binlow = 10.0;
  double binwidth = 10. ;
  for(int i = 0;i < hpteff->GetNbinsX(); i++)
      {
        if(pt > (binlow + i * binwidth) && pt < (binlow + (i+1)*binwidth))
       theeff = hpteff->GetBinContent(i+1);
      }
  return theeff;
}                                                                                                                       



// ------------ method called once each job just before starting event loop  ------------
void 
DimuonsZselAnalyzer::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DimuonsZselAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DimuonsZselAnalyzer);
