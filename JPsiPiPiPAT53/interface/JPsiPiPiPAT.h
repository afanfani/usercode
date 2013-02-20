// -*- C++ -*-
//
// Package:    JPsiPiPiPAT
// Class:      JPsiPiPiPAT
// 
/**\class JPsiPiPiPAT JPsiPiPiPAT.cc myAnalyzers/JPsiPiPiPAT/src/JPsiPiPiPAT.cc

 Description: <one line class summary>
Make rootTuple for JPsiPiPi reconstruction

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//

#ifndef _JPsiPiPiPAT_h
#define _JPsiPiPiPAT_h

// system include files
#include <memory>
#include <map>
#include <string>

// user include files
#include "../interface/VertexReProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoVertex/V0Producer/interface/V0Producer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
//
// class decleration
//

using std::vector;
using namespace edm;
using namespace reco;
using namespace std;

class JPsiPiPiPAT : public edm::EDAnalyzer {
public:
  explicit JPsiPiPiPAT(const edm::ParameterSet&);
  ~JPsiPiPiPAT();
  
private:
  
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

	InvariantMassFromVertex massCalculator;
	
	const reco::DeDxDataValueMap * energyLoss;
	//const reco::TrackRef & trackref;
	Int_t iexception_dedx;

	virtual double
	getSigmaOfLogdEdx(double logde);
	virtual float
	getEnergyLoss(const reco::TrackRef & track);
	virtual double
	nsigmaofdedx(const reco::TrackRef & track, double &theo, double &sigma);
	virtual double
	getLogdEdx(double bg);
	virtual double
	GetMass(const reco::TrackRef & track);
    	virtual double
    	GetMassMC(const reco::TrackRef & track);

    bool isAbHadron(int pdgID);
    bool isAMixedbHadron(int pdgID, int momPdgID);
    std::pair<int, float> findCandMCInfo(reco::GenParticleRef genCand);

  
  // ----------member data ---------------------------
  std::string proccessName_;
  HLTConfigProvider hltConfig_;

  edm::InputTag hlTriggerResults_;
  std::map<std::string,int> *HLTTrig; // HLT trigger prescale for accepted paths

  edm::InputTag inputGEN_;
  std::string vtxSample;
  bool doMC;
  int MCParticle;
  bool doJPsiMassCost;
  bool skipPsi2S;
  bool samesign;

  int MuPixHits_c;
  int MuSiHits_c;
  double MuNormChi_c;
  double MuD0_c;

  double JMaxM_c;
  double JMinM_c;
  double  PsiMaxM_c;
  double PsiMinM_c;
  int PiSiHits_c;
  double PiPt_c;
  double JPiPiDR_c;
  double XPiPiDR_c;
  bool UseXDr_c;	
  double JPiPiMax_c;
  double JPiPiMin_c;
	
  bool resolveAmbiguity_; 
  bool addXlessPrimaryVertex_;
  vector<string>      TriggersForMatching_;
  vector<string>      FiltersForMatching_;
  int  MatchingTriggerResult[50];
  bool Debug_;
  double Chi_Track_;
  
  std::string DeDxEstimator_;    
  bool pileupInfoPresent_;

//AF test
  int nevt_PreSelNo;
  int nevt_PreSelMuon;
  int nevt_PreSelJPsi;

  TTree* X_One_Tree_;

  	unsigned int        runNum, evtNum, lumiNum,NumberVtx;
  	vector<unsigned int>* trigRes;
  	vector<std::string>* trigNames;
  	//  vector<unsigned int>* trigPreScl;
  	vector<unsigned int>* L1TT;
  	vector<std::string>* MatchTriggerNames;

  	unsigned int	nX, nJ, nMu, nMC, nMCAll, nB;

  	float	priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxChiNorm, priVtxChi, priVtxCL; // primary vertex in the event

        // X candidate :
  	vector<float>	*xMass, *xVtxCL, *xVtxC2, *xPx, *xPy, *xPz;
  	vector<double>	*xPxE, *xPyE, *xPzE,*xMass_err,*NumberTrPV, *WeightOfVertex, *SumOfPtPV;
  	vector<float>	*xDecayVtxX, *xDecayVtxY, *xDecayVtxZ;      // 4-tracks vertex
  	vector<double>	*xDecayVtxXE, *xDecayVtxYE, *xDecayVtxZE;
  	vector<float>	*PriVtxXCorrX, *PriVtxXCorrY, *PriVtxXCorrZ;     
  	vector<double> 	*PriVtxXCorrEX, *PriVtxXCorrEY, *PriVtxXCorrEZ;
  	vector<float>	*PriVtxXCorrC2, *PriVtxXCorrCL;
  	vector<double> 	*xLxyPV, *xCosAlpha, *xCTauPV;  // Lifetime wrt PV without X (removing muons and pions) and closest in z 
        vector<double>  *xLxyBS, *xCosAlphaBS, *xCTauBS; 
  	vector<double>	*xCTauPVE, *xCTauBSE;
  	vector<double>	*xLxyPVX, *xCosAlphaX, *xCTauPVX, *xCTauPVEX; //  Lifetime wrt PV with smaller longitudinal X impact parameter
  	vector<float>	*xCTauPVX_3D, *xCTauPVX_3D_err;
  
  	vector<int>	*JIndex; // JPsi index
  	vector<int>	*pipIdx, *pimIdx; // pion index
        vector<int>     *mupIdx, *mumIdx; // muon index
        // JPsi :
	vector<float>	*JMass, *JVtxCL, *JVtxC2, *JPx, *JPy, *JPz;
	vector<float>	*JDecayVtxX, *JDecayVtxY, *JDecayVtxZ;
	vector<float>	*JDecayVtxXE, *JDecayVtxYE, *JDecayVtxZE;
	vector<bool>	*JPsiMuonTrigMatch;
        // negative muons (m) and positive muon (p):
	vector<float>	*mumPx, *mumPy, *mumPz, *mumfChi2;
	vector<float>	*mupPx, *mupPy, *mupPz, *mupfChi2;
	vector<int>	*mumfNDF, *mupfNDF, *jtype;
        vector<float>   *muPx, *muPy, *muPz, *muChi2, *muGlChi2,  *mufHits;
    	vector<unsigned int>	*muKey;
	vector<bool>    *muFirstBarrel, *muFirstEndCap;
	vector<float>	*muD0E, *muDzVtxErr,*muD0, *muDz, *muDzVtx, *muDxyVtx,*muGlDzVtx, *muGlDxyVtx;
  	vector<int>     *muNDF, *muPhits, *muShits, *muLayersTr, *muLayersPix, *muGlNDF, *muGlMuHits, *muGlMatchedStation,*muQual,*muType;
  	vector<int>     *muTrack;
  	vector<float>	*muCharge;
        vector<float>   *MC_muPx, *MC_muPy, *MC_muPz, *MC_muCharge, *MC_muPdgId;
        vector<int>     *MC_muMother1, *MC_muMother2;


    	vector<float>  *BMass, *BPx,  *BPy, *BPz, *BPxE, *BPyE, *BPzE, *BVtxCL, *BVtxC2, *BDecayVtxX, *BDecayVtxY, *BDecayVtxZ, *BDecayVtxXE, *BDecayVtxYE, *BDecayVtxZE;
    	vector<int>    *XIndex_forB, *KaonIdx;
  //vector<float>       *mcPx, *mcPy, *mcPz, *mcE; 
  // 0: X or Psi(2S)
  // 1: J/Psi
  // 2&3 mu+ and mu-
  // 4&5 pi+ and pi-
    	vector<float>   *MCPx, *MCPy, *MCPz;
    	vector<int>     *MCPdgIdAll;
    	vector<float>   *MCPhi,*MCPhiAll,*MC2piPhi, *MCeta, *MCetaAll, *MCE, *MCPt, *MCP, *MCVx, *MCVy, *MCVz, *MCtheta; 
    	vector<int>     *MCPdgId,  *MCnDecay, *MCNDaughters, *MCParent, *MCStatus, *MC_Real_Mom;
    	vector <float>  *ppdlTrue;

        // tracks:
  	vector<float>   *trPx, *trPy, *trPz, *trE, *trPt, *trPtErr;
        vector<float>   *trD0, *trD0E, *trCharge;
        vector<int>     *trQualityHighPurity,*trQualityTight;
	vector<int>     *trNDF, *trPhits, *trShits, *trLayersTr, *trLayersPix;
	vector<float>   *trChi2;
	vector<float>   *trDzVtx, *trDxyVtx;
	vector<double>* tr_nsigdedx;
    	vector<float>* 	tr_dedx, *tr_dedxMass, *tr_theo, *tr_sigma, *tr_dedxMassMC;
        
  	vector<float>   *MC_trPx, *MC_trPy, *MC_trPz,*MC_trphi, *MC_treta, *MC_trPdgId ,*MC_trE;
  	vector<int>	*MC_trMother1, *MC_trMother2;
    
        // muon and pions momenta after kinematic fit
	vector<float>       *fpi1Px, *fpi1Py, *fpi1Pz, *fpi1E;
	vector<float>       *fpi2Px, *fpi2Py, *fpi2Pz, *fpi2E;
	vector<float>       *fmu1Px, *fmu1Py, *fmu1Pz, *fmu1E;
	vector<float>       *fmu2Px, *fmu2Py, *fmu2Pz, *fmu2E;
};

#endif
