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


// system include files
#include <memory>

// user include files
#include "../interface/JPsiPiPiPAT.h"
#include "../interface/VertexReProducer.h"
//#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

//for 53x 
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "TMath.h"
#include "Math/VectorUtil.h"


//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//
JPsiPiPiPAT::JPsiPiPiPAT(const edm::ParameterSet& iConfig)
:
hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN",edm::InputTag("genParticles"))),
vtxSample(iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices"))), 
doMC( iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", false) ),
MCParticle( iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443) ), //20443 X, 100443 Psi(2S), 9120443 X from B
doJPsiMassCost( iConfig.getUntrackedParameter<bool>("DoJPsiMassConstraint", false) ),

skipPsi2S( iConfig.getUntrackedParameter<bool>("SkipPsi2S", true) ),
samesign( iConfig.getUntrackedParameter<bool>("SameSign", false) ),

MuPixHits_c(iConfig.getUntrackedParameter<int>("MinNumMuPixHits", 0)),
MuSiHits_c(iConfig.getUntrackedParameter<int>("MinNumMuSiHits", 0)),
MuNormChi_c(iConfig.getUntrackedParameter<double>("MaxMuNormChi2", 1000)),
MuD0_c(iConfig.getUntrackedParameter<double>("MaxMuD0", 1000)),

JMaxM_c(iConfig.getUntrackedParameter<double>("MaxJPsiMass", 3.2)),
JMinM_c(iConfig.getUntrackedParameter<double>("MinJPsiMass", 2.7)),
PsiMaxM_c(iConfig.getUntrackedParameter<double>("MaxPsi2SMass", 3.7)),
PsiMinM_c(iConfig.getUntrackedParameter<double>("MinPsi2SMass", 3.5)),
PiSiHits_c(iConfig.getUntrackedParameter<int>("MinNumTrSiHits", 0)),
PiPt_c(iConfig.getUntrackedParameter<double>("MinTrPt", 0)),
JPiPiDR_c(iConfig.getUntrackedParameter<double>("JPsiPiPiMaxDR", 1)),
XPiPiDR_c(iConfig.getUntrackedParameter<double>("XCandPiPiMaxDR", 1.1)),
UseXDr_c(iConfig.getUntrackedParameter<bool>("UseXDr", false)),

JPiPiMax_c(iConfig.getUntrackedParameter<double>("JPsiPiPiMaxMass", 50)),
JPiPiMin_c(iConfig.getUntrackedParameter<double>("JPsiPiPiMinMass", 0)),
resolveAmbiguity_(iConfig.getUntrackedParameter<bool>("resolvePileUpAmbiguity",true)),
addXlessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addXlessPrimaryVertex",true)),
TriggersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
FiltersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching")),
Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",false)),
Chi_Track_(iConfig.getUntrackedParameter<double>("Chi2NDF_Track", 10)),
DeDxEstimator_(iConfig.getUntrackedParameter<std::string>("DeDxEstimator",std::string("dedxHarmonic2"))),
pileupInfoPresent_(iConfig.getUntrackedParameter<bool>("pileupInfoPresent",true)),

//AF test
nevt_PreSelNo(0),
nevt_PreSelMuon(0),
nevt_PreSelJPsi(0),


X_One_Tree_(0),runNum(0), evtNum(0), lumiNum(0),NumberVtx(0),

trigRes(0), trigNames(0), L1TT(0),MatchTriggerNames(0), //trigPreScl(0),

nX(0), nJ(0), nMu(0), nMC(0),nMCAll(0),
priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxChiNorm(0), priVtxChi(0),priVtxCL(0),

xMass(0), xVtxCL(0), xVtxC2(0), xPx(0), xPy(0), xPz(0), xPxE(0), xPyE(0), xPzE(0), xMass_err(0),NumberTrPV(0), WeightOfVertex(0), SumOfPtPV(0),
xDecayVtxX(0), xDecayVtxY(0), xDecayVtxZ(0), xDecayVtxXE(0), xDecayVtxYE(0), xDecayVtxZE(0),
PriVtxXCorrX(0), PriVtxXCorrY(0), PriVtxXCorrZ(0), PriVtxXCorrEX(0), PriVtxXCorrEY(0), PriVtxXCorrEZ(0),  PriVtxXCorrC2(0),PriVtxXCorrCL(0),
xLxyPV(0), xCosAlpha(0), xCTauPV(0),  xLxyBS(0), xCosAlphaBS(0), xCTauBS(0), xCTauPVE(0), xCTauBSE(0),
xLxyPVX(0), xCosAlphaX(0), xCTauPVX(0), xCTauPVEX(0),xCTauPVX_3D(0),xCTauPVX_3D_err(0),
JIndex(0), pipIdx(0), pimIdx(0),// mupxIdx(0), mumxIdx(0),

JMass(0), JVtxCL(0), JVtxC2(0), JPx(0), JPy(0), JPz(0),
JDecayVtxX(0), JDecayVtxY(0), JDecayVtxZ(0), JDecayVtxXE(0), JDecayVtxYE(0), JDecayVtxZE(0),
mupIdx(0), mumIdx(0),
JPsiMuonTrigMatch(0),

mumPx(0), mumPy(0), mumPz(0), mumfChi2(0),mupPx(0), mupPy(0), mupPz(0), mupfChi2(0),mumfNDF(0), mupfNDF(0), jtype(0),

muPx(0), muPy(0), muPz(0), muD0(0), muDz(0), muChi2(0), muGlChi2(0),
mufHits(0), muFirstBarrel(0), muFirstEndCap(0),
muDzVtx(0), muDxyVtx(0),
MC_muPx(0), MC_muPy(0), MC_muPz(0), MC_muCharge(0), MC_muPdgId(0),MC_muMother1(0),MC_muMother2(0),
muNDF(0), muGlNDF(0),muPhits(0), muShits(0),muLayersTr(0),muLayersPix(0),muD0E(0), muDzVtxErr(0),muKey(0), muGlMuHits(0), muGlMatchedStation(0), muGlDzVtx(0), muGlDxyVtx(0),
muType(0), muQual(0),
 
muTrack(0), muCharge(0),
BMass(0),BPx(0),BPy(0),BPz(0),BPxE(0),BPyE(0),BPzE(0),BVtxCL(0),BVtxC2(0),BDecayVtxX(0),BDecayVtxY(0),BDecayVtxZ(0),BDecayVtxXE(0),BDecayVtxYE(0),BDecayVtxZE(0),XIndex_forB(0),KaonIdx(0),nB(0),

//  mcPx(0), mcPy(0), mcPz(0), mcE(0),
MCPx(0), MCPy(0), MCPz(0),  MCPdgIdAll(0), MCPhi(0),MCPhiAll(0),  MC2piPhi(0), MCeta(0),  MCetaAll(0),MCE(0), MCPt(0), MCP(0), MCVx(0), MCVy(), MCVz(), MCtheta(),   MCPdgId(0),MCnDecay(0), MCNDaughters(0), MCParent(0),MCStatus(0),MC_Real_Mom(0),ppdlTrue(0),

trPx(0), trPy(0), trPz(0), trE(0),
trPt(0), trPtErr(0),
trNDF(0), trPhits(0), trShits(0), trLayersTr(0),trLayersPix(0),trChi2(0),
trDzVtx(0), trDxyVtx(0),


MC_trPx(0), MC_trPy(0), MC_trPz(0),MC_trphi(0), MC_treta(0), MC_trPdgId(0), MC_trE(0),MC_trMother1(0),MC_trMother2(0),
trD0(0), trD0E(0), trCharge(0), tr_nsigdedx(0), tr_dedx(0), tr_dedxMass(0), tr_theo(0), tr_sigma(0), tr_dedxMassMC(0),
trQualityHighPurity(0), trQualityTight(0),

//AF

fpi1Px(0), fpi1Py(0), fpi1Pz(0), fpi1E(0),
fpi2Px(0), fpi2Py(0), fpi2Pz(0), fpi2E(0),
fmu1Px(0), fmu1Py(0), fmu1Pz(0), fmu1E(0),
fmu2Px(0), fmu2Py(0), fmu2Pz(0), fmu2E(0)
{
    //string DoubleMu0 = "HLT_DoubleMu0"; 
    //TriggersForMatching.push_back(DoubleMu0);
    //string DoubleMu0_Quark = "HLT_DoubleMuO_Quarkonium_V1";
    //TriggersForMatching.push_back(DoubleMu0_Quark);
    //now do what ever initialization is needed
}


JPsiPiPiPAT::~JPsiPiPiPAT()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
    if (Debug_){
        cout<<"Preseletion No = "<<nevt_PreSelNo<<endl;
        cout<<"Preseletion Muon = "<<nevt_PreSelMuon<<endl;
        cout<<"Preseletion JPsi = "<<nevt_PreSelJPsi<<endl;
    }
	
}


//
// member functions
//

// ------------ method called to for each event  ------------
void JPsiPiPiPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
    // get event content information
	
    runNum = iEvent.id().run();
    evtNum = iEvent.id().event();
    lumiNum = iEvent.id().luminosityBlock();
    bool hasRequestedTrigger = false;
    ESHandle<MagneticField> bFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    
    
    if (Debug_) nevt_PreSelNo++;
    
    // first get HLT results
    map<string,int> HLTPreScaleMap;
    edm::Handle<edm::TriggerResults> hltresults;
    try {
        iEvent.getByLabel(hlTriggerResults_,hltresults);
    }
    catch ( ... ) {
        cout << "Couldn't get handle on HLT Trigger!" << endl;
    }
    if (!hltresults.isValid()) {
        cout << "No Trigger Results!" << endl;
    } 
    else {
        int ntrigs=hltresults->size();
        if (ntrigs==0){
            cout << "No trigger name given in TriggerResults of the input " << endl;
        } 
		
        // get hold of trigger names - based on TriggerResults object!
        edm::TriggerNames triggerNames_;
        triggerNames_ = iEvent.triggerNames(*hltresults);
        
        int ntriggers = TriggersForMatching_.size();
        for (int MatchTrig=0; MatchTrig<ntriggers; MatchTrig++){
            MatchingTriggerResult[MatchTrig]= 0;
        }
        
        for (int itrig=0; itrig< ntrigs; itrig++) {
            string trigName = triggerNames_.triggerName(itrig);
            int hltflag = (*hltresults)[itrig].accept();
            // if (Debug_) cout <<"TRIGGER"<<trigName << " " <<hltflag <<endl;
            trigRes->push_back(hltflag);
            trigNames->push_back(trigName);
            
            int ntriggers = TriggersForMatching_.size();
            for (int MatchTrig=0; MatchTrig<ntriggers; MatchTrig++){
                if (TriggersForMatching_[MatchTrig]== triggerNames_.triggerName(itrig)){
                    MatchingTriggerResult[MatchTrig]= hltflag;
                    if (hltflag==1) hasRequestedTrigger =true;
                    break;
                }
            }
        }
        //int ntriggers = TriggersForMatching_.size();
        for (int MatchTrig=0; MatchTrig<ntriggers; MatchTrig++){
            //cout << TriggersForMatching_[MatchTrig]<<endl;
            MatchTriggerNames->push_back(TriggersForMatching_[MatchTrig]);
        }
        
        //
        // Get HLT map : triggername associated with its prescale, saved only for accepted trigger
        //
        for (unsigned int i=0; i<triggerNames_.size(); i++){
            if ( hltresults->accept(i) ) { //  save trigger info only for accepted paths
                // get the prescale from the HLTConfiguration, initialized at beginRun
                int prescale = hltConfig_.prescaleValue(iEvent,iSetup,triggerNames_.triggerNames().at(i));
                //std::cout<<" HLT===> "<<triggerNames_.triggerNames().at(i)<<" prescale ="<<prescale<<std::endl;
                HLTPreScaleMap[triggerNames_.triggerNames().at(i)] = prescale;
            }
        }
        HLTTrig = &HLTPreScaleMap; // store in the branch
        
    } // end valid trigger
    
    // get L1 trigger info
	

    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
	
    edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
    const DecisionWord dWord = gtRecord->decisionWord();  
    
    const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
    for(unsigned int l1i=0; l1i!=ttWord.size(); ++l1i){
        L1TT->push_back(ttWord.at(l1i));
    }

    Vertex thePrimaryV;
    Vertex theBeamSpotV;
    math::XYZPoint RefVtx;
	
    reco::BeamSpot beamSpot;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
    if ( beamSpotHandle.isValid() ) {beamSpot = *beamSpotHandle; 
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
    else cout << "No beam spot available from EventSetup" << endl;
	
	
    Handle<VertexCollection> recVtxs;
    iEvent.getByLabel(vtxSample, recVtxs);
    unsigned int nVtxTrks = 0;
    NumberVtx= recVtxs->size();
    if (Debug_){
        cout<<"Vertexis:"<< NumberVtx <<endl;
    }

    if ( recVtxs->begin() != recVtxs->end() ) {
        if (addXlessPrimaryVertex_ || resolveAmbiguity_){
            thePrimaryV = Vertex(*(recVtxs->begin())); // default primary vertex, should be based on sum of tracks's pT
        }
        else {
            for (reco::VertexCollection::const_iterator vtx = recVtxs->begin();
                 vtx != recVtxs->end(); ++vtx){
                if (nVtxTrks < vtx->tracksSize() ) {
                    nVtxTrks = vtx->tracksSize();
                    thePrimaryV = Vertex(*vtx); // based on max number of tracks
                }
				
            }
        }
    }
    else {
        thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
	
    edm::ESHandle<TransientTrackBuilder> theTTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
    KalmanVertexFitter vtxFitter(true);
    
    RefVtx = thePrimaryV.position();   //reference primary vertex choosen
	
    priVtxX=(thePrimaryV.position().x());
    priVtxY=(thePrimaryV.position().y());
    priVtxZ=(thePrimaryV.position().z());
    priVtxXE=(thePrimaryV.xError());	
    priVtxYE=(thePrimaryV.yError());
    priVtxZE=(thePrimaryV.zError());
    priVtxChiNorm=(thePrimaryV.normalizedChi2());
    priVtxChi=thePrimaryV.chi2();
    priVtxCL = ChiSquaredProbability((double)(thePrimaryV.chi2()),(double)(thePrimaryV.ndof()));
	
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // try reconstruction without fitting modules
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
    Handle< vector<pat::GenericParticle> > thePATTrackHandle;
    iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);
	
    Handle< vector<pat::Muon> > thePATMuonHandle;
    //iEvent.getByLabel("cleanPatMuons", thePATMuonHandle);
    iEvent.getByLabel("patMuonsWithTrigger", thePATMuonHandle);
	
	Handle < reco::DeDxDataValueMap > elossCollection;
	energyLoss = 0;
	iexception_dedx = 0;
	try
	{
		//  iEvent.getByLabel("energyLoss", "energyLossStrHits", elossCollection);
		//  iEvent.getByLabel("dedxMedian", elossCollection);
		iEvent.getByLabel(DeDxEstimator_, elossCollection);
		energyLoss = elossCollection.product();
	} catch (cms::Exception& ex)
	{
		if (evtNum < 100)
			edm::LogError("Analyzer")
					<< "Warning can't get collection with label : elossCollection";
		iexception_dedx = 1;
	}

	
    if (Debug_){
        cout << "starting event with " << thePATTrackHandle->size() << " tracks, and " << thePATMuonHandle->size() << " muons" << endl;
    }
    if (thePATMuonHandle->size()>=2  && hasRequestedTrigger) {
         if (Debug_) cout <<"Accept event with 2 mu and TRIGGER"<<endl;

    // if (thePATMuonHandle->size()>=2){
        //filling track tree
        for ( vector<pat::GenericParticle>::const_iterator iTr = thePATTrackHandle->begin();
             iTr != thePATTrackHandle->end(); ++iTr ) {
            pat::GenericParticle tr = *iTr;
            trPx->push_back(tr.px());
            trPy->push_back(tr.py());
            trPz->push_back(tr.pz());
            trE->push_back(tr.energy());
            trPt->push_back(tr.pt());
            trPtErr->push_back(tr.track()->ptError()); 		

            trNDF->push_back(tr.track()->ndof());
            trPhits->push_back(tr.track()->hitPattern().numberOfValidPixelHits());
            //cout << " Pixel new? " << tr.track()->hitPattern().pixelLayersWithMeasurement() << " pixel old "<<tr.track()->hitPattern().numberOfValidPixelHits()<< " barrel "<< tr.track()->hitPattern().hasValidHitInFirstPixelBarrel() << " endcap " << tr.track()->hitPattern().hasValidHitInFirstPixelEndcap() << endl;
            trShits->push_back(tr.track()->hitPattern().numberOfValidStripHits());
            trLayersTr->push_back(tr.track()->hitPattern().trackerLayersWithMeasurement());
            trLayersPix->push_back(tr.track()->hitPattern().pixelLayersWithMeasurement());
            trChi2->push_back(tr.track()->chi2());
			
            trD0->push_back(tr.track()->d0());
            trD0E->push_back(tr.track()->d0Error());
            trCharge->push_back(tr.charge());
            
            double theo = 0.;
			double sigma = 0.;
			//trackref = tr.track();

			tr_nsigdedx->push_back(nsigmaofdedx(tr.track(), theo, sigma));
			tr_dedx->push_back(getEnergyLoss(tr.track()));
			tr_dedxMass->push_back(GetMass(tr.track()));
            tr_dedxMassMC->push_back(GetMassMC(tr.track()));
			tr_theo->push_back(theo);
			tr_sigma->push_back(sigma);
            

//  Track quality            
//// loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, looseSetWithPV=5, highPuritySetWithPV=6
  bool ishighPurity = tr.track()->quality(reco::TrackBase::highPurity);
  trQualityHighPurity->push_back(ishighPurity);
  trQualityTight->push_back(tr.track()->quality(reco::TrackBase::tight));

//
            float hits = (1.0*tr.track()->found() )/ (tr.track()->found()+ tr.track()->lost() + tr.track()->trackerExpectedHitsInner().numberOfHits() + tr.track()->trackerExpectedHitsOuter().numberOfHits());
            //cout <<"Track "<< hits << " found "<< tr.track()->found() <<" lost " << tr.track()->lost() << " inner "<<tr.track()->trackerExpectedHitsInner().numberOfHits() << " outer "<< tr.track()->trackerExpectedHitsOuter().numberOfHits() << " strip hits "<< tr.track()->hitPattern().numberOfValidStripHits() << "pixel hits" << tr.track()->hitPattern().numberOfValidPixelHits() << " pt "<< tr.px()<<endl;
            trDzVtx->push_back(tr.track()->dz(RefVtx));
            trDxyVtx->push_back(tr.track()->dxy(RefVtx));
   
			
            if (doMC){	 
                if  ((iTr->genParticleRefs().size())>0) {
                    MC_trPx->push_back(iTr->genParticle(0)->px());
                    MC_trPy->push_back(iTr->genParticle(0)->py());
                    MC_trPz->push_back(iTr->genParticle(0)->pz());
                    MC_trphi->push_back(iTr->genParticle(0)->phi());
                    MC_treta->push_back(iTr->genParticle(0)->eta());
                    MC_trE->push_back(iTr->genParticle(0)->energy());
                    MC_trPdgId->push_back(iTr->genParticle(0)->pdgId());
                    if (iTr->genParticle(0)->numberOfMothers()>0){
                        MC_trMother1->push_back(iTr->genParticle(0)->mother()->pdgId());
                        if (iTr->genParticle(0)->mother()->numberOfMothers()>0){
                            MC_trMother2->push_back(iTr->genParticle(0)->mother()->mother()->pdgId());
                        }
                        else {
                            MC_trMother2->push_back(0);
                        }
                    }
                    else {
                        MC_trMother1->push_back(0);
                    }
                    
                }
                
            }
			
        }
		
        //get X and J cands
        for ( std::vector<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin();
             iMuonP != thePATMuonHandle->end(); ++iMuonP) {
            // push back all muon information
            ++nMu;
            const reco::Muon* rmu = dynamic_cast<const reco::Muon * >(iMuonP->originalObject());
            muPx->push_back(rmu->px());
            muPy->push_back(rmu->py());
            muPz->push_back(rmu->pz());
            muCharge->push_back(rmu->charge());
			
            if (doMC){	
                if  ((iMuonP->genParticleRefs().size())>0){
                    MC_muPx->push_back(iMuonP->genParticle(0)->px());
                    MC_muPy->push_back(iMuonP->genParticle(0)->py());
                    MC_muPz->push_back(iMuonP->genParticle(0)->pz());
                    MC_muCharge->push_back(iMuonP->genParticle(0)->charge());
                    MC_muPdgId->push_back(iMuonP->genParticle(0)->pdgId());
                    if (iMuonP->genParticle(0)->numberOfMothers()>0){
                        MC_muMother1->push_back(iMuonP->genParticle(0)->mother()->pdgId());
                        if (iMuonP->genParticle(0)->mother()->numberOfMothers()>0){
                            MC_muMother2->push_back(iMuonP->genParticle(0)->mother()->mother()->pdgId());
                        }
                        else {
                            MC_muMother2->push_back(0);
                        }
                    }
                    else {
                        MC_muMother1->push_back(0);
                    }
                    
                }
            } // end (doMC)
			
			
            if(rmu->track().isNull()){ // rmu->track() returns innerTrack();
                if (Debug_)		 cout << "no track for " << std::distance(thePATMuonHandle->begin(), iMuonP) << " filling defaults" << endl;
                muD0->push_back(0);
                muDz->push_back(0);
                muChi2->push_back(0);
                muNDF->push_back(-1);
                muPhits->push_back(0);
                muShits->push_back(0);
                muLayersTr->push_back(0);
                muLayersPix->push_back(0);
                muDzVtx->push_back(0);
                muDxyVtx->push_back(0);
                mufHits->push_back(0);
                muFirstBarrel->push_back(0);
                muFirstEndCap->push_back(0);
                muD0E->push_back(0);
		muDzVtxErr->push_back(0);
		muKey->push_back(0);
            }
            else {
                muD0->push_back(rmu->track()->d0());
                muDz->push_back(rmu->track()->dz());
                muChi2->push_back(rmu->track()->chi2());
                muNDF->push_back(rmu->track()->ndof());
                muPhits->push_back(rmu->track()->hitPattern().numberOfValidPixelHits());
                muShits->push_back(rmu->track()->hitPattern().numberOfValidStripHits());
                muLayersTr->push_back(rmu->track()->hitPattern().trackerLayersWithMeasurement());
		muLayersPix->push_back(rmu->track()->hitPattern().pixelLayersWithMeasurement());	
                muDzVtx->push_back(rmu->track()->dz(RefVtx));
                muDxyVtx->push_back(rmu->track()->dxy(RefVtx));
                mufHits->push_back((1.0*rmu->track()->found())/ (rmu->track()->found()+ rmu->track()->lost() + rmu->track()->trackerExpectedHitsInner().numberOfHits() + rmu->track()->trackerExpectedHitsOuter().numberOfHits() ) );
                if (Debug_)	  cout <<"mu found " <<rmu->track()->found()<< " " <<(1.0*rmu->track()->found())/ (rmu->track()->found()+ rmu->track()->lost() + rmu->track()->trackerExpectedHitsInner().numberOfHits() + rmu->track()->trackerExpectedHitsOuter().numberOfHits() ) <<endl;
                muFirstBarrel->push_back(rmu->track()->hitPattern().hasValidHitInFirstPixelBarrel());
                muFirstEndCap->push_back(rmu->track()->hitPattern().hasValidHitInFirstPixelEndcap());
                muD0E->push_back(rmu->track()->d0Error());
				muDzVtxErr->push_back(rmu->track()->dzError());
				muKey->push_back(rmu->track().key());
                
            }
            
            if(rmu->globalTrack().isNull()){ // rmu->globalTrack() returns globalTrack();
                if (Debug_)		 cout << "no track for " << std::distance(thePATMuonHandle->begin(), iMuonP) << " filling defaults" << endl;
                muGlChi2->push_back(0);
                muGlNDF->push_back(-1);
                muGlMuHits->push_back(0);
                muGlMatchedStation->push_back(0);
                muGlDxyVtx->push_back(0);
                muGlDzVtx->push_back(0);
            }
            else {
                muGlChi2->push_back(rmu->globalTrack()->chi2());
                muGlNDF->push_back(rmu->globalTrack()->ndof());
                muGlMuHits->push_back(rmu->globalTrack()->hitPattern().numberOfValidMuonHits());
                muGlMatchedStation->push_back(rmu->numberOfMatchedStations());
                muGlDxyVtx->push_back(rmu->muonBestTrack()->dxy(RefVtx));
                muGlDzVtx->push_back(rmu->muonBestTrack()->dz(RefVtx));
            }
            muType->push_back(rmu->type());
            int qm = 0;
            for(int qi=1; qi!= 24; ++qi){
                if(muon::isGoodMuon(*rmu, muon::SelectionType(qi))){
                    qm += 1<<qi;
                }
            }
            muQual->push_back(qm);
            muTrack->push_back(-1);// not implemented yet
			
            //check for mu+
            //       if (iMuonP->charge() != 1) continue;
            TrackRef muTrackP = iMuonP->track();
            if ( muTrackP.isNull() ) {
                if (Debug_)	 cout << "continue due to no track ref" << endl;
                continue;
            }
            if (rmu->track()->hitPattern().numberOfValidPixelHits() < MuPixHits_c 
                || rmu->track()->hitPattern().numberOfValidStripHits() < MuSiHits_c
                || rmu->track()->chi2()/rmu->track()->ndof() > MuNormChi_c
                || fabs(rmu->track()->dxy(RefVtx)) > MuD0_c){
                continue;
            }
            if (Debug_) cout <<"end mu +" <<endl; 
            //next check for mu-
            for ( std::vector<pat::Muon>::const_iterator iMuonM = iMuonP+1;
                 iMuonM != thePATMuonHandle->end(); ++iMuonM) {
                if (Debug_) cout <<"start mu -" <<endl;
                
                //	 if (iMuonM->charge() != -1) continue;
                if(iMuonM->charge() * iMuonP->charge() > 0) continue;
                TrackRef muTrackM = iMuonM->track();
                if ( muTrackM.isNull() ) {
                    //	   cout << "continue from no track ref" << endl;
                    continue;
                }
                const reco::Muon* rmu2 = dynamic_cast<const reco::Muon * >(iMuonM->originalObject());	 
                if(muon::overlap(*rmu, *rmu2)){
                    if (Debug_)   cout << "continued from overlapped muons" << endl;
                    continue;
                }
                if (rmu2->track()->hitPattern().numberOfValidPixelHits() < MuPixHits_c 
                    || rmu2->track()->hitPattern().numberOfValidStripHits() < MuSiHits_c
                    || rmu2->track()->chi2()/rmu->track()->ndof() > MuNormChi_c
                    || fabs(rmu2->track()->dxy(RefVtx)) > MuD0_c){
                    continue;
                }
                if (Debug_) cout <<"bothmuok" <<endl;
                if (Debug_) nevt_PreSelMuon++;
                
                //Get The J/Psi information				
                TransientTrack muonPTT(muTrackP, &(*bFieldHandle) );
                TransientTrack muonMTT(muTrackM, &(*bFieldHandle) );
				
                KinematicParticleFactoryFromTransientTrack pFactory;
				
                //The mass of a muon and the insignificant mass sigma 
                //to avoid singularities in the covariance matrix.
                ParticleMass muon_mass = 0.10565837; //pdg mass
                float muon_sigma = muon_mass*1.e-6;
				
                //initial chi2 and ndf before kinematic fits.
                float chi = 0.;
                float ndf = 0.;
                vector<RefCountedKinematicParticle> muonParticles;
                muonParticles.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
                muonParticles.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
				
                KinematicParticleVertexFitter fitter;   
                RefCountedKinematicTree psiVertexFitTree;
                psiVertexFitTree = fitter.fit(muonParticles); 
                if (!psiVertexFitTree->isValid()) {
                    if (Debug_)	   std::cout << "caught an exception in the psi vertex fit" << std::endl;
                    continue; 
                }
                psiVertexFitTree->movePointerToTheTop();
				
                RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
                RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
				
                psiVertexFitTree->movePointerToTheFirstChild();
                RefCountedKinematicParticle muPCandMC = psiVertexFitTree->currentParticle();
                psiVertexFitTree->movePointerToTheNextChild();
                RefCountedKinematicParticle muMCandMC = psiVertexFitTree->currentParticle();
				
                KinematicParameters psiMupKP = muPCandMC->currentState().kinematicParameters();
                KinematicParameters psiMumKP = muMCandMC->currentState().kinematicParameters();
				
                //Fill the jtree vectors
                JMass->push_back( psi_vFit_noMC->currentState().mass() );
				
                JDecayVtxX->push_back( psi_vFit_vertex_noMC->position().x() );
                JDecayVtxY->push_back( psi_vFit_vertex_noMC->position().y() );
                JDecayVtxZ->push_back( psi_vFit_vertex_noMC->position().z() );
				
                JDecayVtxXE->push_back( sqrt(psi_vFit_vertex_noMC->error().cxx()) );
                JDecayVtxYE->push_back( sqrt(psi_vFit_vertex_noMC->error().cyy()) );
                JDecayVtxZE->push_back( sqrt(psi_vFit_vertex_noMC->error().czz()) );
                JVtxCL->push_back( ChiSquaredProbability((double)(psi_vFit_vertex_noMC->chiSquared()),(double)(psi_vFit_vertex_noMC->degreesOfFreedom())) );
                JVtxC2->push_back( psi_vFit_vertex_noMC->chiSquared() );
				
                JPx->push_back( psiMumKP.momentum().x() + psiMupKP.momentum().x() );
                JPy->push_back( psiMumKP.momentum().y() + psiMupKP.momentum().y() );
                JPz->push_back( psiMumKP.momentum().z() + psiMupKP.momentum().z() );
                mupIdx->push_back(std::distance(thePATMuonHandle->begin(), iMuonP)); 
                mumIdx->push_back(std::distance(thePATMuonHandle->begin(), iMuonM));
				
                mumPx->push_back(psiMumKP.momentum().x());
                mumPy->push_back(psiMumKP.momentum().y());
                mumPz->push_back(psiMumKP.momentum().z());
                mupPx->push_back(psiMupKP.momentum().x());
                mupPy->push_back(psiMupKP.momentum().y());
                mupPz->push_back(psiMupKP.momentum().z());
				
                mumfChi2->push_back(muMCandMC->chiSquared());
                mumfNDF->push_back(muMCandMC->degreesOfFreedom());
                mupfChi2->push_back(muPCandMC->chiSquared());
                mupfNDF->push_back(muPCandMC->degreesOfFreedom());
                
                if (Debug_) cout<<"fine J"<<endl;
                int  dimuonType = 0;   //0 nothing,  1 J/Psi  , 2 Psi(2S)   
                if  ( psi_vFit_noMC->currentState().mass() > JMinM_c && psi_vFit_noMC->currentState().mass() < JMaxM_c)     {
                    dimuonType =1;
                }
                
                if  ( psi_vFit_noMC->currentState().mass() > PsiMinM_c && psi_vFit_noMC->currentState().mass() < PsiMaxM_c)     {
                    dimuonType =2;
                }
                //cout << dimuonType <<endl;
                jtype->push_back(dimuonType);
                
                
                
                int ntriggers =  TriggersForMatching_.size();
                for (int MatchTrig=0; MatchTrig < ntriggers; MatchTrig++){
                    
                    //cout << "result " << MatchingTriggerResult[MatchTrig] << " " << TriggersForMatching_[MatchTrig] <<endl;
                    if (MatchingTriggerResult[MatchTrig]!=0){
                        const pat::TriggerObjectStandAloneCollection mu1HLTMatches = iMuonP->triggerObjectMatchesByFilter(FiltersForMatching_[MatchTrig]  );
                        const pat::TriggerObjectStandAloneCollection mu2HLTMatches = iMuonM->triggerObjectMatchesByFilter(FiltersForMatching_[MatchTrig]  );
                        bool pass1 = mu1HLTMatches.size() > 0;
                        bool pass2 = mu2HLTMatches.size() > 0;
                        
                        //cout << "Uno " << pass1 << " " << TriggersForMatching_[MatchTrig] <<endl;
                        //cout << "Due " << pass2 << " " << TriggersForMatching_[MatchTrig] <<endl;
                        if ((pass1) && (pass2))
                            JPsiMuonTrigMatch->push_back(true);
                        
                        else
                            JPsiMuonTrigMatch->push_back(false);
                    }
                    else
                        JPsiMuonTrigMatch->push_back(false);
                }
				
                ++nJ;
                muonParticles.clear();
				
                if ( thePATTrackHandle->size() < 2) {
                    continue;
                }
                               
                if (dimuonType ==0){
                    continue;
                }

                if ( skipPsi2S && (dimuonType ==2) ){ // skip Psi(2S)pipi decay
                    continue;
                }                

                if (Debug_) nevt_PreSelJPsi++;
                
                //next check tracks for pi+
                //bool match = false;
                //	 /*	 Comment out all of the X candidate stuff
                for ( vector<pat::GenericParticle>::const_iterator iTrackP = thePATTrackHandle->begin();
                     iTrackP != thePATTrackHandle->end(); ++iTrackP ) {
					
                    //check track doesn't overlap with the Psi candidate tracks
                    if (iTrackP->track().key()==rmu->track().key() || iTrackP->track().key()==rmu2->track().key())
                    {//cout << "My Match"<<endl;
                        continue;}
					
                    if( ( iTrackP->track()->chi2()/ iTrackP->track()->ndof() > Chi_Track_) ||  iTrackP->pt() < PiPt_c){
                        continue;
                    }
					
                    //next check tracks for pi-
                    for ( vector<pat::GenericParticle>::const_iterator iTrackM = iTrackP+1;//thePATTrackHandle->begin();
                         iTrackM != thePATTrackHandle->end(); ++iTrackM ) {
                        // Select only tracks with opposite charge      
                        if( !samesign && (iTrackP->charge() * iTrackM->charge() > 0) ) continue;
                        // Select only tracks with same charge :     
                        if( samesign && (iTrackP->charge() * iTrackM->charge() < 0) ) continue;
                        //check track doesn't overlap with the Psi candidate tracks
                        if (iTrackM->track().key()==rmu->track().key() || iTrackM->track().key()==rmu2->track().key())
                        {//cout << "My Match"<<endl;
                            continue;}
						
                        //if( iTrackM->track()->hitPattern().numberOfValidStripHits() < PiSiHits_c
                        if ((iTrackM->track()->chi2()/ iTrackM->track()->ndof()  > Chi_Track_)
                            || iTrackM->pt() < PiPt_c){
                            continue;
                        }
						
                        math::XYZTLorentzVector jp = (rmu->p4() + rmu2->p4());
                        math::XYZTLorentzVector Xp = (rmu->p4() + rmu2->p4() + iTrackP->p4() + iTrackM->p4()); 
                        float jppDR = sqrt(pow(jp.eta()-iTrackP->p4().eta(),2) + pow(jp.phi()-iTrackP->p4().phi(),2));
                        float jpmDR = sqrt(pow(jp.eta()-iTrackM->p4().eta(),2) + pow(jp.phi()-iTrackM->p4().phi(),2));
                        
                        float XppDR = sqrt(pow(Xp.eta()-iTrackP->p4().eta(),2) + pow(Xp.phi()-iTrackP->p4().phi(),2));
                        float XpmDR = sqrt(pow(Xp.eta()-iTrackM->p4().eta(),2) + pow(Xp.phi()-iTrackM->p4().phi(),2));
                        
                        if (UseXDr_c){
                            if (XppDR > XPiPiDR_c || XpmDR > XPiPiDR_c ) continue;
                        }
                        else {
                            if(jppDR > JPiPiDR_c || jpmDR > JPiPiDR_c ) continue;
                        }
                        
                        
                        if((iTrackP->p4() + iTrackM->p4() + jp).M() >= JPiPiMax_c || (iTrackP->p4() + iTrackM->p4() + jp).M() < JPiPiMin_c) continue;
                        
                        
                        //have two oppositely charged muons, and two oppositely charged pions try to vertex them
                        ///Memmory Problem is somewhere between here and a point idicated below	     
                        TransientTrack pionPTT(iTrackP->track(), &(*bFieldHandle) );
                        TransientTrack pionMTT(iTrackM->track(), &(*bFieldHandle) );
						
                        const ParticleMass pion_mass = 0.13957018;
                        float pion_sigma = pion_mass*1.e-6;
						
                        // Do mass constraint for JPsi cand and do mass constrained vertex fit
                        // creating the constraint with a small sigma to put in the resulting covariance 
                        // matrix in order to avoid singularities
						
                        vector<RefCountedKinematicParticle> vFitMCParticles;
                        vFitMCParticles.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
                        vFitMCParticles.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
                        vFitMCParticles.push_back(pFactory.particle(pionPTT,pion_mass,chi,ndf,pion_sigma));
                        vFitMCParticles.push_back(pFactory.particle(pionMTT,pion_mass,chi,ndf,pion_sigma));
                        
						
                        RefCountedKinematicParticle xCandMC;
                        RefCountedKinematicVertex xDecayVertexMC;
                        RefCountedKinematicTree vertexFitTree;
						
                        if(doJPsiMassCost){  // dimuon mass constraint
                            //	       cout << "Doing mass constraint fit" << endl;
                            ParticleMass psi_mass = 3.096916;
                            ParticleMass psi2_mass = 3.68609;
                            MultiTrackKinematicConstraint *  j_psi_c;
                            if (dimuonType ==1){ // constrain to JPsi mass
                                j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
                            }
                            if (dimuonType ==2){ // constrain to Psi(2S) mass
                                j_psi_c = new  TwoTrackMassKinematicConstraint(psi2_mass);
                            }
                            
                            KinematicConstrainedVertexFitter kcvFitter;
                            vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
                            if (!vertexFitTree->isValid()) {
                                //	       std::cout << "caught an exception in the X vertex fit with MC" << std::endl;
                                continue;
                            }
							
                            vertexFitTree->movePointerToTheTop();
                            xCandMC = vertexFitTree->currentParticle();
                            xDecayVertexMC = vertexFitTree->currentDecayVertex();
                        } else {
                            KinematicParticleVertexFitter kcvFitter;
                            vertexFitTree = kcvFitter.fit(vFitMCParticles);
                            if (!vertexFitTree->isValid()) {
                                // 		 std::cout << "caught an exception in the X vertex fit with MC" << std::endl;
                                continue;
                            }
							
                            vertexFitTree->movePointerToTheTop();
                            xCandMC = vertexFitTree->currentParticle();
                            xDecayVertexMC = vertexFitTree->currentDecayVertex();
                        }
                        if (!xDecayVertexMC->vertexIsValid()){
                            //	       cout << "X MC fit vertex is not valid" << endl;
                            continue;
                        }
                        
                        double xVtxProb = ChiSquaredProbability((double)(xDecayVertexMC->chiSquared()),(double)(xDecayVertexMC->degreesOfFreedom()));
                        if ( xVtxProb<0.005 ) continue;
                        
                        
                        if ( xCandMC->currentState().mass() > 20 ) {continue;}
                        
                        if (Debug_){
                            cout<<"Starting mass error"<<endl;
                        }
                        vector<double> PionMasses;
                        PionMasses.push_back( 0.13957018 );
                        PionMasses.push_back( 0.13957018 );
                        vector<TransientTrack> t_tksPions;
                        t_tksPions.push_back(theTTBuilder->build(*iTrackP->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
                        t_tksPions.push_back(theTTBuilder->build(*iTrackM->track()));
                        
                        //TransientVertex myPionVertex = vtxFitter.vertex(t_tksPions);
                        
                        CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tksPions );
                        Measurement1D MassWErr = massCalculator.invariantMass( VtxForInvMass, PionMasses );
                        
                        xMass_err->push_back(MassWErr.error());
                        if (Debug_){
                            cout<<"Ending mass error"<<endl;
                        }
                        
                        //Memory problem is above this line
						
                        //	     /* Comment out fill of vectors

                        // fill X() candidate variables now
                        xMass->push_back(xCandMC->currentState().mass());
                        xPx->push_back(xCandMC->currentState().globalMomentum().x());
                        xPy->push_back(xCandMC->currentState().globalMomentum().y());
                        xPz->push_back(xCandMC->currentState().globalMomentum().z());
						
                        xPxE->push_back( sqrt( xCandMC->currentState().kinematicParametersError().matrix()(3,3) ) );
                        xPyE->push_back( sqrt( xCandMC->currentState().kinematicParametersError().matrix()(4,4) ) );
                        xPzE->push_back( sqrt( xCandMC->currentState().kinematicParametersError().matrix()(5,5) ) );
						
                        xVtxCL->push_back( ChiSquaredProbability((double)(xDecayVertexMC->chiSquared()),(double)(xDecayVertexMC->degreesOfFreedom())) );
                        xVtxC2->push_back( xDecayVertexMC->chiSquared() );
                        xDecayVtxX->push_back((*xDecayVertexMC).position().x());
                        xDecayVtxY->push_back((*xDecayVertexMC).position().y());
                        xDecayVtxZ->push_back((*xDecayVertexMC).position().z());
                        xDecayVtxXE->push_back(sqrt((*xDecayVertexMC).error().cxx()));
                        xDecayVtxYE->push_back(sqrt((*xDecayVertexMC).error().cyy()));
                        xDecayVtxZE->push_back(sqrt((*xDecayVertexMC).error().czz()));
                        vertexFitTree->movePointerToTheFirstChild();
                        RefCountedKinematicParticle mu1 =  vertexFitTree->currentParticle();
                        vertexFitTree->movePointerToTheNextChild();
                        RefCountedKinematicParticle mu2 = vertexFitTree->currentParticle();
                        vertexFitTree->movePointerToTheNextChild();
                        RefCountedKinematicParticle pi1 = vertexFitTree->currentParticle();
                        vertexFitTree->movePointerToTheNextChild();
                        RefCountedKinematicParticle pi2 = vertexFitTree->currentParticle();
                        // fill pions after fit
                        fpi1Px->push_back(pi1->currentState().globalMomentum().x());
                        fpi1Py->push_back(pi1->currentState().globalMomentum().y());
                        fpi1Pz->push_back(pi1->currentState().globalMomentum().z());
                        fpi1E->push_back(pi1->currentState().kinematicParameters().energy());
                        fpi2Px->push_back(pi2->currentState().globalMomentum().x());
                        fpi2Py->push_back(pi2->currentState().globalMomentum().y());
                        fpi2Pz->push_back(pi2->currentState().globalMomentum().z());
                        fpi2E->push_back(pi2->currentState().kinematicParameters().energy());
						
                        // fill muons after fit
                        fmu1Px->push_back(mu1->currentState().globalMomentum().x());
                        fmu1Py->push_back(mu1->currentState().globalMomentum().y());
                        fmu1Pz->push_back(mu1->currentState().globalMomentum().z());
                        fmu1E->push_back(mu1->currentState().kinematicParameters().energy());
                        fmu2Px->push_back(mu2->currentState().globalMomentum().x());
                        fmu2Py->push_back(mu2->currentState().globalMomentum().y());
                        fmu2Pz->push_back(mu2->currentState().globalMomentum().z());
                        fmu2E->push_back(mu2->currentState().kinematicParameters().energy());
						
		        Vertex theOriginalPV = thePrimaryV;
					
                        vector<TransientVertex> pvs;
                        if (addXlessPrimaryVertex_) {
                            VertexReProducer revertex(recVtxs, iEvent);
                            Handle<TrackCollection> pvtracks;   
                            iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
                            Handle<BeamSpot>        pvbeamspot;
                            iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);
                            if (pvbeamspot.id() != beamSpotHandle.id()) edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
                            const reco::Muon *rmu_1 = dynamic_cast<const reco::Muon *>(iMuonP->originalObject());
                            const reco::Muon *rmu_2 = dynamic_cast<const reco::Muon *>(iMuonM->originalObject());
                            if (rmu_1 != 0 && rmu_2 != 0 && rmu_1->track().id() == pvtracks.id() && rmu_2->track().id() == pvtracks.id() && iTrackP->track().id() ==  pvtracks.id() && iTrackM->track().id() ==  pvtracks.id()) { 
                                TrackCollection XLess;
                                XLess.reserve(pvtracks->size());
                                for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
                                    if (i == rmu_1->track().key()) continue;
                                    if (i == rmu_2->track().key()) continue;
                                    if (i == iTrackP->track().key()) continue;
                                    if (i == iTrackM->track().key()) continue;
                                    XLess.push_back((*pvtracks)[i]);
                                }
                                pvs = revertex.makeVertices(XLess, *pvbeamspot, iSetup) ;
                                if (!pvs.empty()) {
                                    Vertex XLessPV = Vertex(pvs.front()); // get the first one
                                    thePrimaryV = XLessPV;
                                }
                            }
                        }
						
                        Vertex theOtherV = thePrimaryV; 
                        //cout<<" choose PV ="<<endl;
                        if (resolveAmbiguity_) {
                            float minDz = 999999.;
                            float minDzTrack = 999999.;
                            if (!addXlessPrimaryVertex_) {
                                for(VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv){
                                    float deltaZ = fabs((*xDecayVertexMC).position().z() - itv->position().z()) ;
                                    if ( deltaZ < minDz ) {
                                        minDz = deltaZ;    
                                        thePrimaryV = Vertex(*itv);
                                    }
                                }
                            } else {
                                for(vector<TransientVertex>::iterator itv2 = pvs.begin(), itvend2 = pvs.end(); itv2 != itvend2; ++itv2){
                                    float deltaZ = fabs((*xDecayVertexMC).position().z() - itv2->position().z()) ;
                                    math::XYZPoint CurrentVtx = Vertex(*itv2).position();
                                    float DzTrack = fabs(xCandMC->refittedTransientTrack().track().dz(CurrentVtx));
                                    //cout<<" z(X) - z(vtx) ="<<deltaZ<<endl;
                                    //cout<<" X dz ="<<DzTrack<<endl; 
                                    if ( deltaZ < minDz ) {
                                        minDz = deltaZ;    
                                        Vertex XLessPV = Vertex(*itv2); 
                                        theOtherV = XLessPV;
                                        //cout<<" z(X) - z(vtx) min="<<minDz<<endl; 
                                    }
                                    if ( DzTrack<minDzTrack){
                                        minDzTrack = DzTrack;
                                        Vertex XLessPV = Vertex(*itv2);
                                        theOtherV = XLessPV;
                                        //cout<<" X dz min="<<minDzTrack<<endl;
                                    }

                                }
                            }
                        } 
                        
                        
                        Vertex TheOtherVertex3D = thePrimaryV; 
                        //cout<<" choose PV ="<<endl;
                        if (resolveAmbiguity_) {
                            float minDz = 999999.;
                            float minDzTrack = 999999.;
                            if (!addXlessPrimaryVertex_) {
                                for(VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv){
                                    float deltaZ = fabs((*xDecayVertexMC).position().z() - itv->position().z()) ;
                                    if ( deltaZ < minDz ) {
                                        minDz = deltaZ;    
                                        TheOtherVertex3D = Vertex(*itv);
                                    }
                                }
                            } else {
                                for(vector<TransientVertex>::iterator itv2 = pvs.begin(), itvend2 = pvs.end(); itv2 != itvend2; ++itv2){
                                     VertexDistance3D a3d;
  				     float deltaZ   = a3d.distance(Vertex(*itv2), Vertex(*xDecayVertexMC)).value();
                                     if ( deltaZ < minDz ) {
                                        minDz = deltaZ;    
                                        Vertex XLessPV = Vertex(*itv2); 
                                        TheOtherVertex3D = XLessPV;
                                        //cout<<" z(X) - z(vtx) min="<<minDz<<endl; 
                                    }

                                }
                            }
                        } 
                        
						
                        PriVtxXCorrX->push_back(thePrimaryV.position().x());
                        PriVtxXCorrY->push_back(thePrimaryV.position().y());
                        PriVtxXCorrZ->push_back(thePrimaryV.position().z());
                        PriVtxXCorrEX->push_back(thePrimaryV.xError());
                        PriVtxXCorrEY->push_back(thePrimaryV.yError());
                        PriVtxXCorrEZ->push_back(thePrimaryV.zError());
                        PriVtxXCorrCL->push_back( ChiSquaredProbability((double)(thePrimaryV.chi2()),(double)(thePrimaryV.ndof())) );
                        PriVtxXCorrC2->push_back(thePrimaryV.chi2());
						
                        //Correction tu muDz e muDxy with new PV
                        
                        //muDzVtx[(std::distance(thePATMuonHandle->begin(), iMuonP))] = rmu->track()->dz(thePrimaryV));
                        //muDxyVtx[(std::distance(thePATMuonHandle->begin(), iMuonP))] = rmu->track()->dxy(thePrimaryV));
                        //muDzVtx[(std::distance(thePATMuonHandle->begin(), iMuonM))] = rmu2->track()->dz(thePrimaryV));
                        //muDxyVtx[(std::distance(thePATMuonHandle->begin(), iMuonM))] = rmu2->track()->dxy(thePrimaryV));
                        
		       	float Pt_Tot =0;
		       	int N_tracks = 0;
                        try{				
				//cout<<"Counting Tracks"<<endl;
				for(reco::Vertex::trackRef_iterator itVtx = theOriginalPV.tracks_begin(); itVtx != theOriginalPV.tracks_end(); itVtx++){
				   const reco::Track& track1 = **itVtx;
				   if(!track1.quality(reco::TrackBase::highPurity)) continue;
				   if (track1.pt() < 0.6)  continue;
				   if ((track1.ptError()/ track1.pt())>0.1) continue;
     
                                   float Tr_X_DR = sqrt(pow(track1.eta()-xCandMC->currentState().globalMomentum().eta(),2) + pow(track1.phi()-xCandMC->currentState().globalMomentum().phi(),2));
                                   if (Tr_X_DR<0.7){
                                    Pt_Tot+=track1.pt();
                                    N_tracks++;
                                   }
                            
				}
			} catch (std::exception & err) {//std::cout << "counting tracks from PVertex failed " << std::endl; 
						 return ; }
			if (Debug_) cout <<" counting tracks from PVertex (with hardcoded track selection) "<< Pt_Tot << " " << N_tracks<<endl;
                        NumberTrPV->push_back(N_tracks);
                        WeightOfVertex->push_back( 0);
                        SumOfPtPV->push_back( Pt_Tot);

                        
                        TVector3 vtx;  //4-track fit vertex
                        vtx.SetXYZ((*xDecayVertexMC).position().x(), (*xDecayVertexMC).position().y(), 0);
						
                        TVector3 pvtx; 
                        pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
						
                        VertexDistanceXY vdistXY;
                        
                        TVector3 pperp(xCandMC->currentState().globalMomentum().x(),xCandMC->currentState().globalMomentum().y(), 0);
                        
                        AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
						
                        TVector3 vdiff = vtx - pvtx;
                        double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                        Measurement1D distXY = vdistXY.distance(Vertex(*xDecayVertexMC), Vertex(thePrimaryV));
                        double ctauPV = distXY.value()*cosAlpha*xCandMC->currentState().mass()/pperp.Perp();
						
						
                        GlobalError v1e = (Vertex(*xDecayVertexMC)).error();
                        GlobalError v2e = thePrimaryV.error();
                        AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
                       // double ctauErrPV = sqrt(vXYe.similarity(vpperp))*xCandMC->currentState().mass()/(pperp.Perp2());
                        double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*xCandMC->currentState().mass()/(pperp.Perp2());
                        float lxyPV = vdiff.Dot(pperp)/pperp.Mag();
						
                        xCosAlpha->push_back(cosAlpha);
                        xCTauPV->push_back(ctauPV);
                        xCTauPVE->push_back(ctauErrPV);
                        xLxyPV->push_back(lxyPV);
			
                        // Lifetime wrt PV with smaller longitudinal X impact parameter 
                        pvtx.SetXYZ(theOtherV.position().x(), theOtherV.position().y(), 0);
                        vdiff = vtx - pvtx;
                        cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                        distXY = vdistXY.distance(Vertex(*xDecayVertexMC), Vertex(theOtherV));
                        double ctauPVX = distXY.value()*cosAlpha*xCandMC->currentState().mass()/pperp.Perp();
                        GlobalError v1eX = (Vertex(*xDecayVertexMC)).error();
                        GlobalError v2eX = theOtherV.error();
                        AlgebraicSymMatrix33 vXYeX = v1eX.matrix()+ v2eX.matrix();
                        //double ctauErrPVX = sqrt(vXYeX.similarity(vpperp))*xCandMC->currentState().mass()/(pperp.Perp2());
                        double ctauErrPVX = sqrt(ROOT::Math::Similarity(vpperp,vXYeX))*xCandMC->currentState().mass()/(pperp.Perp2());
                        float lxyPVX = vdiff.Dot(pperp)/pperp.Mag();

                        xCosAlphaX->push_back(cosAlpha);
                        xCTauPVX->push_back(ctauPVX);
                        xCTauPVEX->push_back(ctauErrPVX);
                        xLxyPVX->push_back(lxyPVX);
                        
                        
                        
                        
                        
                        VertexDistance3D a3d;
  						float Dist3DPV     = a3d.distance(TheOtherVertex3D, Vertex(*xDecayVertexMC)).value();
  						float Dist3DPV_err    = a3d.distance(TheOtherVertex3D, Vertex(*xDecayVertexMC)).error();
  						xCTauPVX_3D->push_back(Dist3DPV);
  						xCTauPVX_3D_err->push_back(Dist3DPV_err);
						//cout << Dist3DPV << " " << Dist3DPV_err << endl; 
			
                        //Lifetime BS
                        pvtx.SetXYZ(theBeamSpotV.position().x(), theBeamSpotV.position().y(), 0);
                        vdiff = vtx - pvtx;
                        cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                        distXY = vdistXY.distance(Vertex(*xDecayVertexMC), Vertex(theBeamSpotV));
                        double ctauBS = distXY.value()*cosAlpha*xCandMC->currentState().mass()/pperp.Perp();
                        GlobalError v1eB = (Vertex(*xDecayVertexMC)).error();
                        GlobalError v2eB = theBeamSpotV.error();
                        AlgebraicSymMatrix33 vXYeB = v1eB.matrix()+ v2eB.matrix();
                        //double ctauErrBS = sqrt(vXYeB.similarity(vpperp))*xCandMC->currentState().mass()/(pperp.Perp2());
                        double ctauErrBS = sqrt(ROOT::Math::Similarity(vpperp,vXYeB))*xCandMC->currentState().mass()/(pperp.Perp2());
                        float lxyBS = vdiff.Dot(pperp)/pperp.Mag();
						
                        xCosAlphaBS->push_back(cosAlpha);
                        xCTauBS->push_back(ctauBS);
                        xCTauBSE->push_back(ctauErrBS);
                        xLxyBS->push_back(lxyBS);
                        
						
                        JIndex->push_back(nJ-1);
                        pipIdx->push_back(std::distance(thePATTrackHandle->begin(), iTrackP));
                        pimIdx->push_back(std::distance(thePATTrackHandle->begin(), iTrackM));
                        //	     */  // Comment out Fill of X vectors
						
                        nX++;

                        // Create a state J/Psi + pi + pi + K
                        //cout << xCandMC->currentState().mass() <<endl;
                        //cout <<ctauPV*10<< endl;
//                        if  ( (fabs(xCandMC->currentState().mass()-3.871)<0.03) && (ctauPV*10>0.08) && (dimuonType ==1)){ 
//                      restrict to X() region and to decay with JPsi
                          if  ( (xCandMC->currentState().mass())>3.8 && (xCandMC->currentState().mass())<3.95 && (dimuonType ==1)){
                            for ( vector<pat::GenericParticle>::const_iterator iTrackK = thePATTrackHandle->begin();
                                 iTrackK != thePATTrackHandle->end(); ++iTrackK ) {
                                if (iTrackK->track().key()==rmu->track().key() || iTrackK->track().key()==rmu2->track().key() ||  (iTrackK==iTrackM)  ||  (iTrackK==iTrackP))
                                {//cout << "My Match"<<endl;
                                    continue;}
                                
                                if( ( iTrackK->track()->chi2()/ iTrackK->track()->ndof() > Chi_Track_)  || iTrackK->pt() < PiPt_c){
                                    //cout << "Cinematic"<<endl;
                                    continue;
                                }
                                // restrict to B mass region
                                if(  (iTrackP->p4() + iTrackM->p4() + jp + iTrackK->p4()).M() >= 7.0 || (iTrackP->p4() + iTrackM->p4() + jp + iTrackK->p4()).M() < 4.5) continue;
                                
                                TransientTrack Kaon(iTrackK->track(), &(*bFieldHandle) );
                                
                                vector<RefCountedKinematicParticle> vFitBMeson;
                                const ParticleMass kaon_mass = 0.493677;
                                float kaon_sigma = 0.000016;

                                vFitBMeson.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
                                vFitBMeson.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
                                vFitBMeson.push_back(pFactory.particle(pionPTT,pion_mass,chi,ndf,pion_sigma));
                                vFitBMeson.push_back(pFactory.particle(pionMTT,pion_mass,chi,ndf,pion_sigma));
                                vFitBMeson.push_back(pFactory.particle(Kaon,kaon_mass,chi,ndf,kaon_sigma));
                                
                                RefCountedKinematicParticle BCand;
                                RefCountedKinematicVertex BDecayVertex;
                                RefCountedKinematicTree vertexFitTreeB;
                                
                                //	       cout << "Doing mass constraint fit" << endl;
                                // ParticleMass psi2Smass = 3.68609;
                                
                                 ParticleMass psi_mass_inB = 3.096916;
                                 MultiTrackKinematicConstraint *  conforB = new  TwoTrackMassKinematicConstraint(psi_mass_inB);
                                 KinematicConstrainedVertexFitter fitterB;
                                 vertexFitTreeB = fitterB.fit(vFitBMeson, conforB);
                                
                                //KinematicParticleVertexFitter fitterB;   
                                //RefCountedKinematicTree psiVertexFitTree;
                                //vertexFitTreeB = fitterB.fit(vFitBMeson);
                               // vertexFitTreeB = fitterB.fit(vFitBMeson);
                                
                                if (!vertexFitTreeB->isValid()) {
                                    //	       std::cout << "caught an exception in the X vertex fit with MC" << std::endl;
                                    continue;
                                }

                                vertexFitTreeB->movePointerToTheTop();
                                BCand = vertexFitTreeB->currentParticle();
                                BDecayVertex = vertexFitTreeB->currentDecayVertex();

                                
                                if (!BDecayVertex->vertexIsValid()){
                                    //	       cout << "X MC fit vertex is not valid" << endl;
                                    continue;
                                }

                                if ((ChiSquaredProbability((double)(BDecayVertex->chiSquared()),(double)(BDecayVertex->degreesOfFreedom())))<0.01) {
                                    //cout << " failed chi2 cut in MC fit with chi2 = " << BDecayVertex->chiSquared() << endl;
                                    continue;
                                }
                                
                                if ( BCand->currentState().mass() > 7 ) {continue;}
                                
                                
                                BMass->push_back(BCand->currentState().mass());
                                BPx->push_back(BCand->currentState().globalMomentum().x());
                                BPy->push_back(BCand->currentState().globalMomentum().y());
                                BPz->push_back(BCand->currentState().globalMomentum().z());
                                
                                BPxE->push_back( sqrt( BCand->currentState().kinematicParametersError().matrix()(3,3) ) );
                                BPyE->push_back( sqrt( BCand->currentState().kinematicParametersError().matrix()(4,4) ) );
                                BPzE->push_back( sqrt( BCand->currentState().kinematicParametersError().matrix()(5,5) ) );
                                
                                
                                BVtxCL->push_back( ChiSquaredProbability((double)(BDecayVertex->chiSquared()),(double)(BDecayVertex->degreesOfFreedom())) );
                                BVtxC2->push_back( BDecayVertex->chiSquared() );
                                BDecayVtxX->push_back((*BDecayVertex).position().x());
                                BDecayVtxY->push_back((*BDecayVertex).position().y());
                                BDecayVtxZ->push_back((*BDecayVertex).position().z());
                                BDecayVtxXE->push_back(sqrt((*BDecayVertex).error().cxx()));
                                BDecayVtxYE->push_back(sqrt((*BDecayVertex).error().cyy()));
                                BDecayVtxZE->push_back(sqrt((*BDecayVertex).error().czz()));
                                
                                XIndex_forB->push_back(nX-1);
                                KaonIdx->push_back(std::distance(thePATTrackHandle->begin(), iTrackK));
                                
                                /////
                                
                                nB++;
                                vFitBMeson.clear();
                                
                            }
                            
                           // vFitBMeson.clear();
                            
                        }
                        
                        vFitMCParticles.clear();
                        
                        
                        
                    }// 2nd loop over track (look for pi-)
                }// 1st loop over track (look for pi+)
                //	 */ // End Comment out all of X candidate stuff
            }// 2nd loop over muons (look for mu-)
        }//first loop over muons (look for mu+)
    }//if two muons
	
    //   if (doMC){
    //  mcPx->resize(6, 0);
    //  mcPy->resize(6, 0);
    //   mcPz->resize(6, 0);
    //  mcE->resize(6, 0);
	
    //check variables
	
    if (doMC){
        double genPhi;
        Handle<GenParticleCollection> genParticles;
        iEvent.getByLabel(inputGEN_, genParticles);
        GenParticleCollection::const_iterator genObj = genParticles->begin();
		
        for(; genObj != genParticles->end(); genObj++) {
            int nx3872 = 0;
            MCPdgIdAll->push_back(genObj->pdgId());
            MCPhiAll->push_back(genObj->phi());
            MCetaAll->push_back(genObj->eta());
            nMCAll++;
            
            if (fabs(genObj->pdgId()) == 521){
                //cout <<"qui"<<endl;
                candidate::const_iterator daughterItB  = genObj->begin();
                candidate::const_iterator daughterEndB = genObj->end();
                bool isX =false;
                bool isK = false;
                for(; daughterItB != daughterEndB; ++daughterItB) {
                    if( fabs(daughterItB->pdgId())==MCParticle) {  
                        isX =true;
                    }
                    if( fabs(daughterItB->pdgId())==321) { 
                        isK =true;
                    }
                }
                if (isX && isK){
                    genPhi = (genObj->phi() > 0) ? genObj->phi() : genObj->phi() + 2*M_PI;
                    MCPdgId->push_back(genObj->pdgId());
                    MCPx->push_back(genObj->px());
                    MCPy->push_back(genObj->py());
                    MCPz->push_back(genObj->pz());
                    MCPt->push_back(genObj->pt());
                    MCPhi->push_back(genObj->phi());
                    MC2piPhi->push_back(genPhi);
                    MCeta->push_back(genObj->eta());
                    MCE->push_back(genObj->energy());
                    MCParent->push_back(genObj->mother()->pdgId());
                    MCNDaughters->push_back(genObj->numberOfDaughters());
                    MCnDecay->push_back(1);
                    MCVx->push_back(genObj->vx());
                    MCVy->push_back(genObj->vy());
                    MCVz->push_back(genObj->vz());
                    MCP->push_back(genObj->p());
                    MCtheta->push_back(genObj->theta());
                    MCStatus->push_back(genObj->status());
                    nMC++;
                    
                    candidate::const_iterator daughterItB2  = genObj->begin();
                    candidate::const_iterator daughterEndB2 = genObj->end();
                    for(; daughterItB2 != daughterEndB2; ++daughterItB2) {
                        if ( fabs(daughterItB2->pdgId())==321){
                            genPhi = (daughterItB2->phi() > 0) ? daughterItB2->phi() : daughterItB2->phi() + 2*M_PI;
                            MCPdgId->push_back(daughterItB2->pdgId());
                            MCPx->push_back(daughterItB2->px());
                            MCPy->push_back(daughterItB2->py());
                            MCPz->push_back(daughterItB2->pz());
                            MCPt->push_back(daughterItB2->pt());
                            MCPhi->push_back(daughterItB2->phi());
                            MC2piPhi->push_back(genPhi);
                            MCeta->push_back(daughterItB2->eta());
                            MCE->push_back(daughterItB2->energy());
                            MCParent->push_back(daughterItB2->mother()->pdgId());
                            MCNDaughters->push_back(daughterItB2->numberOfDaughters());
                            MCnDecay->push_back(1);
                            MCVx->push_back(daughterItB2->vx());
                            MCVy->push_back(daughterItB2->vy());
                            MCVz->push_back(daughterItB2->vz());
                            MCP->push_back(daughterItB2->p());
                            MCtheta->push_back(daughterItB2->theta());
                            MCStatus->push_back(daughterItB2->status());
                            nMC++;

                        }
                    }
                    
                }
            }
            
            if( genObj->pdgId() == MCParticle ){
                
                
                if  (genObj->mother()->pdgId() != MCParticle){
                    reco::GenParticleRef test = genObj->daughterRef(0);
                    std::pair<int, float> MCinfo = findCandMCInfo(test->motherRef());
                    //cout<<MCinfo.first<<endl;
                    MC_Real_Mom->push_back(MCinfo.first);
                    ppdlTrue->push_back(MCinfo.second);
                }
                
                // if(( genObj->pdgId() == 20443 )||( genObj->pdgId() == 100443)) {  //chi_c1->X3872
                if (genObj->numberOfDaughters()>1){
                    nx3872++;
                    genPhi = (genObj->phi() > 0) ? genObj->phi() : genObj->phi() + 2*M_PI;
                    MCPdgId->push_back(genObj->pdgId());
                    MCPx->push_back(genObj->px());
                    MCPy->push_back(genObj->py());
                    MCPz->push_back(genObj->pz());
                    MCPt->push_back(genObj->pt());
                    MCPhi->push_back(genObj->phi());
                    MC2piPhi->push_back(genPhi);
                    MCeta->push_back(genObj->eta());
                    MCE->push_back(genObj->energy());
                    MCParent->push_back(genObj->mother()->pdgId());
                    MCNDaughters->push_back(genObj->numberOfDaughters());
                    MCnDecay->push_back(nx3872);
                    MCVx->push_back(genObj->vx());
                    MCVy->push_back(genObj->vy());
                    MCVz->push_back(genObj->vz());
                    MCP->push_back(genObj->p());
                    MCtheta->push_back(genObj->theta());
                    MCStatus->push_back(genObj->status());
                    nMC++;
                    int nJPsi=0;
                    int nJPsid=0;
                    int nrho=0;
                    int nrhod=0;
                    // 	   int pip=0;
                    // 	   int pim=0;
                    // 	   int mup=0;
                    // 	   int mum=0;
                    //cout << "Jpsi"<<endl;
                    candidate::const_iterator daughterIt  = genObj->begin();
                    candidate::const_iterator daughterEnd = genObj->end();
                    for(; daughterIt != daughterEnd; ++daughterIt) {
                        if( fabs(daughterIt->pdgId())==443 ) {  
                            nJPsi++;
                            nJPsid = nJPsid + daughterIt->numberOfDaughters();
                            candidate::const_iterator dIt  = daughterIt->begin();
                            candidate::const_iterator dEnd = daughterIt->end();
                            for(; dIt != dEnd; ++dIt) {
                                // 	     if (dIt->pdgId()==13){
                                // 	       mup++;
                                // 	     }
                                // 	     if (dIt->pdgId()==-13){
                                // 	       mum++;
                                // 	     }
                                // 	     int essow = dIt->pdgId();
                                // 	     if (nJPsid>2)
                                // 	       cout << "componente " << essow<< endl;
                                //dIt->pdgId() <<endl; 
                                genPhi = (dIt->phi() > 0) ? dIt->phi() : dIt->phi() + 2*M_PI;
                                MCPdgId->push_back(dIt->pdgId());
                                MCPx->push_back(dIt->px());
                                MCPy->push_back(dIt->py());
                                MCPz->push_back(dIt->pz());
                                MCPt->push_back(dIt->pt());
								
                                MCVx->push_back(dIt->vx());
                                MCVy->push_back(dIt->vy());
                                MCVz->push_back(dIt->vz());
                                MCP->push_back(dIt->p());
                                MCtheta->push_back(dIt->theta());
                                MCStatus->push_back(dIt->status());
                                MCPhi->push_back(dIt->phi());
                                MC2piPhi->push_back(genPhi);
                                MCeta->push_back(dIt->eta());
                                MCE->push_back(dIt->energy());
                                MCParent->push_back(dIt->mother()->pdgId());
                                MCNDaughters->push_back(dIt->numberOfDaughters());
                                MCnDecay->push_back(nx3872);
                                nMC++;
                            }//loop over J/Psi daughters   
                        }// if jpsi   
                        if( fabs(daughterIt->pdgId())==113 ) {
                            nrho++;
                            nrhod = nrhod + daughterIt->numberOfDaughters();
                            candidate::const_iterator dIt  = daughterIt->begin();
                            candidate::const_iterator dEnd = daughterIt->end();
                            for(; dIt != dEnd; ++dIt) {
                                // 	     if (dIt->pdgId()==211){
                                // 	       pip++;
                                // 	     }
                                // 	     if (dIt->pdgId()==-211){
                                // 	       pim++;
                                // 	     }
                                double genPhi = (dIt->phi() > 0) ? dIt->phi() : dIt->phi() + 2*M_PI;
                                MCPdgId->push_back(dIt->pdgId());
                                MCPx->push_back(dIt->px());
                                MCPy->push_back(dIt->py());
                                MCPz->push_back(dIt->pz());
                                MCPt->push_back(dIt->pt());
                                MCPhi->push_back(dIt->phi());
                                MC2piPhi->push_back(genPhi);
                                MCeta->push_back(dIt->eta());
                                MCE->push_back(dIt->energy());
                                MCParent->push_back(dIt->mother()->pdgId());
                                MCNDaughters->push_back(dIt->numberOfDaughters());
                                MCnDecay->push_back(nx3872);
                                MCVx->push_back(dIt->vx());
                                MCVy->push_back(dIt->vy());
                                MCVz->push_back(dIt->vz());
                                MCP->push_back(dIt->p());
                                MCtheta->push_back(dIt->theta());
                                MCStatus->push_back(dIt->status());
                                nMC++;
                            }// loop over rho daughters 
                        }// if rho   
						
                        double genPhi = (daughterIt->phi() > 0) ? daughterIt->phi() : daughterIt->phi() + 2*M_PI;
                        MCPdgId->push_back(daughterIt->pdgId());
                        MCPx->push_back(daughterIt->px());
                        MCPy->push_back(daughterIt->py());
                        MCPz->push_back(daughterIt->pz());
                        MCPt->push_back(daughterIt->pt());
                        MCPhi->push_back(daughterIt->phi());
                        MC2piPhi->push_back(genPhi);
                        MCeta->push_back(daughterIt->eta());
                        MCE->push_back(daughterIt->energy());
                        MCParent->push_back(daughterIt->mother()->pdgId());
                        MCNDaughters->push_back(daughterIt->numberOfDaughters());
                        MCnDecay->push_back(nx3872);
                        MCVx->push_back(daughterIt->vx());
                        MCVy->push_back(daughterIt->vy());
                        MCVz->push_back(daughterIt->vz());
                        MCP->push_back(daughterIt->p());
                        MCtheta->push_back(daughterIt->theta());
                        MCStatus->push_back(daughterIt->status());
                        nMC++;
                    }//loop over x daughters
                    //cout << "JPsi "<< nJPsi << " d " << nJPsid<< " mum " << mum << " mup "<< mup << " nrho " << nrho<< " pip " <<pip<< " pim " << pim << " nrhpd " << nrhod << endl;
                }// if at least 2 daughters
                //if (mup==1) && (mum==1) && (nJpsi==1) && (nJPsid==2) && nrho==1
            }// if 20443
        }// gen particle loop
    }//if doMC
	
    //fill the tree and clear the vectors
	
    if (nJ > 0 ) {
        //     cout << "Filling Trees" << endl;
        // 		estree_->Fill();
        // 		mutree_->Fill();
        // 		trtree_->Fill();
        // 		xtree_->Fill();
        // 		jtree_->Fill();
        // 		if(doMC) mctree_->Fill();
        X_One_Tree_->Fill();
    }
    if (Debug_){
        cout << "Resetting branches, had " << nX << " X cands  ," << nJ << " J cands and " <<nB <<" B Mother Cand" << endl;
    }
    //    l1_mu3 = 0; l1_2mu3 = 0; l1_muOpen = 0; l1_mu0 = 0; l1_mu7 = 0; l1_mu14 = 0; l1_2muOpen = 0;
    //    hlt_mu3 = 0; hlt_mu5 = 0; hlt_mu7 = 0; hlt_mu9 = 0; hlt_2mu0 = 0; hlt_2mu3 = 0; hlt_2mu3JPsi = 0; hlt_BJPsiMuMu = 0;
    //    hlt_l2mu9 = 0; openhlt_l2mu25 = 0; openhlt_l2mu30 = 0; openhlt_mu11 = 0;
    //    hlt_l1mu14_l1eg10 = 0; hlt_l1mu14_l1jet6u = 0; hlt_l1mu14_l1etm30 = 0;
    //    hlt_mu0tkmu0_jpsi = 0; hlt_mu3tkmu0_jpsi = 0; hlt_mu5tkmu0_jpsi = 0;
    //    hlt_mu0tkmu0_jpsiNC = 0; hlt_mu3tkmu0_jpsiNC = 0; hlt_mu5tkmu0_jpsiNC = 0;
    //    hlt_mu5tk0_jpsi = 0; hlt_mu5tk3_jpsi = 0; hlt_mu5_l2mu0 = 0;
    trigRes->clear(); trigNames->clear(); L1TT->clear(); MatchTriggerNames->clear(); //trigPreScl->clear();MatchTriggerNames
	
    runNum = 0; evtNum=0; lumiNum=0; NumberVtx=0;
    nX = 0; nJ = 0; nMu =0; nMC =0; nB = 0;
	
    priVtxX = 0; priVtxY = 0; priVtxZ = 0; priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxChiNorm = 0; priVtxChi = 0; priVtxCL = 0;
	
    xMass->clear(); xVtxCL->clear(); xVtxC2->clear(); 
    xPx->clear(); xPy->clear(); xPz->clear(); 
    xPxE->clear(); xPyE->clear(); xPzE->clear(); xMass_err->clear(); NumberTrPV->clear(); WeightOfVertex->clear(); SumOfPtPV->clear();
    xDecayVtxX->clear(); xDecayVtxY->clear(); xDecayVtxZ->clear(); 
    xDecayVtxXE->clear(); xDecayVtxYE->clear(); xDecayVtxZE->clear(); 
    PriVtxXCorrX->clear(); PriVtxXCorrY->clear(); PriVtxXCorrZ->clear(); PriVtxXCorrEX->clear(); PriVtxXCorrEY->clear(); PriVtxXCorrEZ->clear();  PriVtxXCorrC2->clear(); PriVtxXCorrCL->clear();
    xLxyPV->clear(); xCosAlpha->clear(); xCTauPV->clear();xLxyBS->clear(); xCosAlphaBS->clear(); xCTauBS->clear(); 
    xCTauPVE->clear(); xCTauBSE->clear();
    xLxyPVX->clear(); xCosAlphaX->clear(); xCTauPVX->clear(); xCTauPVEX->clear();
    xCTauPVX_3D->clear();
  	xCTauPVX_3D_err->clear();

    JIndex->clear();
    pipIdx->clear(); pimIdx->clear();
    // mupxIdx->clear(); mumxIdx->clear();
	
    JMass->clear(); JVtxCL->clear(); JVtxC2->clear(); 
    JPx->clear(); JPy->clear(); JPz->clear();
    JDecayVtxX->clear(); JDecayVtxY->clear(); JDecayVtxZ->clear();
    JDecayVtxXE->clear(); JDecayVtxYE->clear(); JDecayVtxZE->clear();
    mupIdx->clear(); mumIdx->clear();
    JPsiMuonTrigMatch->clear();
	
    mumPx->clear(); mumPy->clear(); mumPz->clear();
    mupPx->clear(); mupPy->clear(); mupPz->clear();
    mumfChi2->clear(); mumfNDF->clear();
    mupfChi2->clear(); mupfNDF->clear(); jtype->clear();
	
    muPx->clear(); muPy->clear(); muPz->clear(); 
    muD0->clear(); muDz->clear(); muChi2->clear(); muGlChi2->clear();
    mufHits->clear(); muFirstBarrel->clear(); muFirstEndCap->clear(); muD0E->clear() ;  muDzVtxErr->clear() ; muKey->clear() ;
    muDzVtx->clear(); muDxyVtx->clear(); muGlDzVtx->clear(); muGlDxyVtx->clear();
    muNDF->clear(); muGlNDF->clear(); muPhits->clear(); muShits->clear(); muLayersTr->clear() ; muLayersPix->clear(); muGlMuHits->clear(); muGlMatchedStation->clear();muType->clear(); 
    muQual->clear(); muTrack->clear(); muCharge->clear();
	
    trPx->clear(); trPy->clear(); trPz->clear(); trE->clear();
    trPt->clear(); trPtErr->clear();
	
    trNDF->clear(); trPhits->clear(); trShits->clear(); trLayersTr->clear(); trLayersPix->clear(); trChi2->clear();
    trD0->clear(); trD0E->clear(); trCharge->clear(); 
    trQualityHighPurity->clear();trQualityTight->clear();
    trDzVtx->clear(); trDxyVtx->clear();
    tr_nsigdedx->clear(); tr_dedx->clear();
    tr_dedxMass->clear(); tr_theo->clear(); tr_sigma->clear(), tr_dedxMassMC->clear();
	//trNSigmaDeDx->clear();
	//trDeDx->clear(); trDeDxError->clear(); trNSM->clear(); trNM->clear();
    
    BMass->clear(); BPx->clear(); BPy->clear(); BPz->clear(); BPxE->clear(); BPyE->clear(); BPzE->clear(); BVtxCL->clear(); BVtxC2->clear(); BDecayVtxX->clear(); BDecayVtxY->clear(); BDecayVtxZ->clear(); BDecayVtxXE->clear(); BDecayVtxYE->clear(); BDecayVtxZE->clear(); XIndex_forB->clear(); KaonIdx->clear();
    
	
    fpi1Px->clear(); fpi1Py->clear(); fpi1Pz->clear(); fpi1E->clear();
    fpi2Px->clear(); fpi2Py->clear(); fpi2Pz->clear(); fpi2E->clear();
    fmu1Px->clear(); fmu1Py->clear(); fmu1Pz->clear(); fmu1E->clear();
    fmu2Px->clear(); fmu2Py->clear(); fmu2Pz->clear(); fmu2E->clear();
	
    if(doMC){
        MC_trPx->clear(); MC_trPy->clear(); MC_trPz->clear(); MC_trphi->clear(); MC_treta->clear(); MC_trE->clear(); MC_trPdgId->clear();MC_trMother1->clear();MC_trMother2->clear();
        MC_muPx->clear(); MC_muPy->clear(); MC_muPz->clear(); MC_muCharge->clear(); MC_muPdgId->clear();MC_muMother1->clear();MC_muMother2->clear();
        MCPx->clear(); MCPy->clear(); MCPz->clear(); MCE->clear(); MCVx->clear(); MCVy->clear(); MCVz->clear(); 
        MCPdgId->clear(); MCPdgIdAll->clear(); MCPhi->clear();  MCPhiAll->clear(); MC2piPhi->clear(); MCeta->clear();  MCetaAll->clear(); MCParent->clear();  
        MCNDaughters->clear();  MCPt->clear();  MCnDecay->clear(); MCtheta->clear(); MCStatus->clear();MC_Real_Mom->clear(); ppdlTrue->clear(); MCP->clear();
    }
}// analyze


// ------------ method called once each job just before starting event loop  ------------
void JPsiPiPiPAT::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup){
    //   bool changed = true;
    //   proccessName_="HLT";
    //   hltConfig_.init(iRun,iSetup,proccessName_,changed);
    //AF:
    bool changed(true);
    if (hltConfig_.init(iRun,iSetup,"HLT",changed)) {
        if (changed) {
            cout<<" HLT config changed"<<endl;
        }
    } else {
        cout<< " HLT config extraction failure with process name HLT ";
    }
    
}


void JPsiPiPiPAT::beginJob()
{
    edm::Service<TFileService> fs;
	
    //estree_ = fs->make<TTree>("eventSummary", "General Event Summary");
    X_One_Tree_ = fs->make<TTree>("X_data", "X(3872) Data");
    
    X_One_Tree_->Branch("TrigRes", &trigRes);
    X_One_Tree_->Branch("TrigNames", &trigNames);
    X_One_Tree_->Branch("MatchTriggerNames", &MatchTriggerNames);
    X_One_Tree_->Branch("HLTTrig","map<string,int>",&HLTTrig);  // map with HLT path and its prescale (stored only for passed trigger)
    //  X_One_Tree__->Branch("TrigPreScl", &trigPreScl);
    X_One_Tree_->Branch("L1TrigRes", &L1TT);
	
    X_One_Tree_->Branch("evtNum",&evtNum,"evtNum/i");
    X_One_Tree_->Branch("runNum",&runNum,"runNum/i");
    X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
    X_One_Tree_->Branch("NumberVtx",&NumberVtx,"NumberVtx/i");
    X_One_Tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
    X_One_Tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
    X_One_Tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
    X_One_Tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
    X_One_Tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
    X_One_Tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
    X_One_Tree_->Branch("priVtxChiNorm",&priVtxChiNorm, "priVtxChiNorm/f");
    X_One_Tree_->Branch("priVtxChi",&priVtxChi, "priVtxChi/f");
    X_One_Tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");
	
    //	X_One_Tree_ = fs->make<TTree>("muons", "Muons");
    X_One_Tree_->Branch("nMu", &nMu, "nMu/i");
    X_One_Tree_->Branch("muPx",&muPx);
    X_One_Tree_->Branch("muPy",&muPy);
    X_One_Tree_->Branch("muPz",&muPz);
    X_One_Tree_->Branch("muD0",&muD0);
    X_One_Tree_->Branch("muDz",&muDz);
    X_One_Tree_->Branch("muChi2",&muChi2);
    X_One_Tree_->Branch("muNDF",&muNDF);
    X_One_Tree_->Branch("muPhits",&muPhits);
    X_One_Tree_->Branch("muShits",&muShits);
    X_One_Tree_->Branch("muLayersTr",&muLayersTr);
    X_One_Tree_->Branch("muLayersPix",&muLayersPix);
    
    X_One_Tree_->Branch("muD0E",&muD0E);
    X_One_Tree_->Branch("muDzVtxErr",&muDzVtxErr);
    X_One_Tree_->Branch("muKey",&muKey);
    
    X_One_Tree_->Branch("muGlChi2",&muGlChi2);
    X_One_Tree_->Branch("muGlNDF",&muGlNDF);
    X_One_Tree_->Branch("muGlMuHits",&muGlMuHits);
    X_One_Tree_->Branch("muGlMatchedStation",&muGlMatchedStation);
    X_One_Tree_->Branch("muGlDzVtx", &muGlDzVtx);
    X_One_Tree_->Branch("muGlDxyVtx", &muGlDxyVtx);
    X_One_Tree_->Branch("muType",&muType);
    X_One_Tree_->Branch("muQual",&muQual);
    X_One_Tree_->Branch("muTrack",&muTrack);
    X_One_Tree_->Branch("muCharge", &muCharge);
   
    X_One_Tree_->Branch("mufHits", &mufHits);
    X_One_Tree_->Branch("muFirstBarrel", &muFirstBarrel);
    X_One_Tree_->Branch("muFirstEndCap", &muFirstEndCap);
    X_One_Tree_->Branch("muDzVtx", &muDzVtx);
    X_One_Tree_->Branch("muDxyVtx", &muDxyVtx);
	
    X_One_Tree_->Branch("TrackPx", &trPx);
    X_One_Tree_->Branch("TrackPy", &trPy);
    X_One_Tree_->Branch("TrackPz", &trPz);
    X_One_Tree_->Branch("TrackEnergy", &trE);
    X_One_Tree_->Branch("TrackPt", &trPt);
    X_One_Tree_->Branch("TrackPtErr", &trPtErr);
    X_One_Tree_->Branch("TrackNDF", &trNDF);
    X_One_Tree_->Branch("TrackPhits", &trPhits);
    X_One_Tree_->Branch("TrackShits", &trShits);
    X_One_Tree_->Branch("TrackLayersTr", &trLayersTr);
    X_One_Tree_->Branch("TrackLayersPix", &trLayersPix);
    X_One_Tree_->Branch("TrackChi2", &trChi2);
    X_One_Tree_->Branch("TrackD0", &trD0);
    X_One_Tree_->Branch("TrackD0Err", &trD0E);
    X_One_Tree_->Branch("TrackCharge", &trCharge);
    X_One_Tree_->Branch("TrackHighPurity", &trQualityHighPurity);
    X_One_Tree_->Branch("TrackTight", &trQualityTight);
    X_One_Tree_->Branch("trDzVtx", &trDzVtx);
    X_One_Tree_->Branch("trDxyVtx", &trDxyVtx);
	
    X_One_Tree_->Branch("tr_nsigdedx", &tr_nsigdedx);
	X_One_Tree_->Branch("tr_dedx", &tr_dedx);
	X_One_Tree_->Branch("tr_dedxMass", &tr_dedxMass);
	X_One_Tree_->Branch("tr_theo", &tr_theo);
	X_One_Tree_->Branch("tr_sigma", &tr_sigma);
    X_One_Tree_->Branch("tr_dedxMassMC", &tr_dedxMassMC);
	
    //	X_One_Tree_ = fs->make<TTree>("Xntuple","Xtomumupipi ntuple");
    X_One_Tree_->Branch("nX",&nX,"nX/i");
    X_One_Tree_->Branch("xMass",&xMass);
    X_One_Tree_->Branch("xVtxCL",&xVtxCL);
    X_One_Tree_->Branch("xVtxC2",&xVtxC2);
    X_One_Tree_->Branch("xPx",&xPx);
    X_One_Tree_->Branch("xPy",&xPy);
    X_One_Tree_->Branch("xPz",&xPz);
    X_One_Tree_->Branch("xPxE",&xPxE);
    X_One_Tree_->Branch("xPyE",&xPyE);
    X_One_Tree_->Branch("xPzE",&xPzE);
    X_One_Tree_->Branch("xMass_err",& xMass_err);
    X_One_Tree_->Branch("NumberTrPV",& NumberTrPV);
    X_One_Tree_->Branch("WeightOfVertex",& WeightOfVertex);
    X_One_Tree_->Branch("SumOfPtPV",& SumOfPtPV);
    
    X_One_Tree_->Branch("xDecayVtxX",&xDecayVtxX);
    X_One_Tree_->Branch("xDecayVtxY",&xDecayVtxY);
    X_One_Tree_->Branch("xDecayVtxZ",&xDecayVtxZ);
    X_One_Tree_->Branch("xDecayVtxXE",&xDecayVtxXE);
    X_One_Tree_->Branch("xDecayVtxYE",&xDecayVtxYE);
    X_One_Tree_->Branch("xDecayVtxZE",&xDecayVtxZE);
	
    X_One_Tree_->Branch("PriVtxXCorrX",&PriVtxXCorrX);
    X_One_Tree_->Branch("PriVtxXCorrY",&PriVtxXCorrY);
    X_One_Tree_->Branch("PriVtxXCorrZ",&PriVtxXCorrZ);
    X_One_Tree_->Branch("PriVtxXCorrEX",&PriVtxXCorrEX);
    X_One_Tree_->Branch("PriVtxXCorrEY",&PriVtxXCorrEY);
    X_One_Tree_->Branch("PriVtxXCorrEZ",&PriVtxXCorrEZ);
    X_One_Tree_->Branch("PriVtxXCorrC2",&PriVtxXCorrC2);
    X_One_Tree_->Branch("PriVtxXCorrCL",&PriVtxXCorrCL);
    X_One_Tree_->Branch("xLxyPV", &xLxyPV);
    X_One_Tree_->Branch("xCosAlpha", &xCosAlpha);
    X_One_Tree_->Branch("xCTauPV", &xCTauPV);
    X_One_Tree_->Branch("xCTauPVE", &xCTauPVE);
    X_One_Tree_->Branch("xLxyBS", &xLxyBS);
    X_One_Tree_->Branch("xCosAlphaBS", &xCosAlphaBS);
    X_One_Tree_->Branch("xCTauBS", &xCTauBS);
    X_One_Tree_->Branch("xCTauBSE", &xCTauBSE);
    X_One_Tree_->Branch("xLxyPVX", &xLxyPVX);
    X_One_Tree_->Branch("xCosAlphaX", &xCosAlphaX);
    X_One_Tree_->Branch("xCTauPVX", &xCTauPVX);
    X_One_Tree_->Branch("xCTauPVEX", &xCTauPVEX);
    X_One_Tree_->Branch("xCTauPVX_3D", &xCTauPVX_3D);
    X_One_Tree_->Branch("xCTauPVX_3D_err", &xCTauPVX_3D_err);
    X_One_Tree_->Branch("JPsiIndex", &JIndex);
    X_One_Tree_->Branch("pipIndex", &pipIdx);
    X_One_Tree_->Branch("pimIndex", &pimIdx);
    X_One_Tree_->Branch("fitpi1Px", &fpi1Px);
    X_One_Tree_->Branch("fitpi1Py", &fpi1Py);
    X_One_Tree_->Branch("fitpi1Pz", &fpi1Pz);
    X_One_Tree_->Branch("fitpi1E", &fpi1E);
    X_One_Tree_->Branch("fitpi2Px", &fpi2Px);
    X_One_Tree_->Branch("fitpi2Py", &fpi2Py);
    X_One_Tree_->Branch("fitpi2Pz", &fpi2Pz);
    X_One_Tree_->Branch("fitpi2E", &fpi2E);
    X_One_Tree_->Branch("fitmu1Px", &fmu1Px);
    X_One_Tree_->Branch("fitmu1Py", &fmu1Py);
    X_One_Tree_->Branch("fitmu1Pz", &fmu1Pz);
    X_One_Tree_->Branch("fitmu1E", &fmu1E);
    X_One_Tree_->Branch("fitmu2Px", &fmu2Px);
    X_One_Tree_->Branch("fitmu2Py", &fmu2Py);
    X_One_Tree_->Branch("fitmu2Pz", &fmu2Pz);
    X_One_Tree_->Branch("fitmu2E", &fmu2E);
	
    //	X_One_Tree_ = fs->make<TTree>("jntuple","Jtomumu ntuple");
	
    X_One_Tree_->Branch("nJ",&nJ,"nJ/i");
    X_One_Tree_->Branch("JMass",&JMass);
    X_One_Tree_->Branch("JVtxCL",&JVtxCL);
    X_One_Tree_->Branch("JVtxC2",&JVtxC2);
    X_One_Tree_->Branch("JPx",&JPx);
    X_One_Tree_->Branch("JPy",&JPy);
    X_One_Tree_->Branch("JPz",&JPz);
    X_One_Tree_->Branch("JDecayVtxX",&JDecayVtxX);
    X_One_Tree_->Branch("JDecayVtxY",&JDecayVtxY);
    X_One_Tree_->Branch("JDecayVtxZ",&JDecayVtxZ);
    X_One_Tree_->Branch("JDecayVtxXE",&JDecayVtxXE);
    X_One_Tree_->Branch("JDecayVtxYE",&JDecayVtxYE);
    X_One_Tree_->Branch("JDecayVtxZE",&JDecayVtxZE);
    X_One_Tree_->Branch("mupIdx", &mupIdx);
    X_One_Tree_->Branch("mumIdx", &mumIdx);
    X_One_Tree_->Branch("mumPx",&mumPx);
    X_One_Tree_->Branch("mumPy",&mumPy);
    X_One_Tree_->Branch("mumPz",&mumPz);
    X_One_Tree_->Branch("mupPx",&mupPx);
    X_One_Tree_->Branch("mupPy",&mupPy);
    X_One_Tree_->Branch("mupPz",&mupPz);
    X_One_Tree_->Branch("mumfChi2",&mumfChi2);
    X_One_Tree_->Branch("mumfNDF",&mumfNDF);
    X_One_Tree_->Branch("mupfChi2",&mupfChi2);
    X_One_Tree_->Branch("mupfNDF",&mupfNDF);
    X_One_Tree_->Branch("jtype",&jtype);
    X_One_Tree_->Branch("JPsiMuonTrigMatch",&JPsiMuonTrigMatch);
    
    
    X_One_Tree_->Branch("BMass",&BMass); 
    X_One_Tree_->Branch("BPx",&BPx);
    X_One_Tree_->Branch("BPy",&BPy);
    X_One_Tree_->Branch("BPz",&BPz);    
    X_One_Tree_->Branch("BPxE",&BPxE);
    X_One_Tree_->Branch("BPyE",&BPyE);
    X_One_Tree_->Branch("BPzE",&BPzE);    
    X_One_Tree_->Branch("BVtxCL",&BVtxCL);
    X_One_Tree_->Branch("BVtxC2",&BVtxC2);
    X_One_Tree_->Branch("BDecayVtxX",&BDecayVtxX);
    X_One_Tree_->Branch("BDecayVtxY",&BDecayVtxY);
    X_One_Tree_->Branch("BDecayVtxZ",&BDecayVtxZ);
    X_One_Tree_->Branch("BDecayVtxXE",&BDecayVtxXE);
    X_One_Tree_->Branch("BDecayVtxYE",&BDecayVtxYE);
    X_One_Tree_->Branch("BDecayVtxZE",&BDecayVtxZE);
    X_One_Tree_->Branch("XIndex_forB",&XIndex_forB); 
    X_One_Tree_->Branch("KaonIdx",&KaonIdx);
    X_One_Tree_->Branch("nB",&nB,"nB/i");
    
    
	
    if(doMC){
        //	X_One_Tree_ = fs->make<TTree>("MCntuple","MC ntuple");
        X_One_Tree_->Branch("nMC",&nMC,"nMC/i");
        X_One_Tree_->Branch("MCPx", &MCPx);
        X_One_Tree_->Branch("MCPy", &MCPy);
        X_One_Tree_->Branch("MCPz", &MCPz);
        X_One_Tree_->Branch("MCEnergy", &MCE);
        X_One_Tree_->Branch("MCEta", &MCeta);
        X_One_Tree_->Branch("MCEtaAll", &MCetaAll);
        X_One_Tree_->Branch("MCPdgIdAll", &MCPdgIdAll);
        X_One_Tree_->Branch("MCPhi", &MCPhi);
        X_One_Tree_->Branch("MCPhiAll", &MCPhiAll);
        X_One_Tree_->Branch("MC2piPhi", &MC2piPhi);
        X_One_Tree_->Branch("MCParent", &MCParent);
        X_One_Tree_->Branch("MCNDaughters", &MCNDaughters);
        X_One_Tree_->Branch("MCPt", &MCPt);
        X_One_Tree_->Branch("MCnDecay", &MCnDecay);
        X_One_Tree_->Branch("MCPdgId", &MCPdgId);
        
		
        X_One_Tree_->Branch("MCVx", &MCVx);
        X_One_Tree_->Branch("MCVy", &MCVy);
        X_One_Tree_->Branch("MCVz", &MCVz);
        X_One_Tree_->Branch("MCMomentum", &MCP);
        X_One_Tree_->Branch("MCTheta", &MCtheta);
        X_One_Tree_->Branch("MCStatus", &MCStatus);
        X_One_Tree_->Branch("MC_Real_Mom", &MC_Real_Mom);
        X_One_Tree_->Branch("ppdlTrue", &ppdlTrue);
        X_One_Tree_->Branch("MCP", &MCP);
        
        
		
        X_One_Tree_->Branch("MC_TrackPx", &MC_trPx);
        X_One_Tree_->Branch("MC_TrackPy", &MC_trPy);
        X_One_Tree_->Branch("MC_TrackPz", &MC_trPz);
        X_One_Tree_->Branch("MC_TrackPhi", &MC_trphi);
        X_One_Tree_->Branch("MC_TrackEta", &MC_treta);
        X_One_Tree_->Branch("MC_TrackPdgId", &MC_trPdgId);
        X_One_Tree_->Branch("MC_TrackEnergy", &MC_trE);
        X_One_Tree_->Branch("MC_TrackMother", &MC_trMother1);
        X_One_Tree_->Branch("MC_TrackMotherMother", &MC_trMother2);
		
		
        X_One_Tree_->Branch("MC_muPx",&MC_muPx);
        X_One_Tree_->Branch("MC_muPy",&MC_muPy);
        X_One_Tree_->Branch("MC_muPz",&MC_muPz);
        X_One_Tree_->Branch("MC_muCharge",&MC_muCharge);
        X_One_Tree_->Branch("MC_muPdgId",&MC_muPdgId);
        X_One_Tree_->Branch("MC_muMother", &MC_muMother1);
        X_One_Tree_->Branch("MC_muMotherMother", &MC_muMother2);
        
		
		
    }
}// begin Job

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiPiPiPAT::endJob() {
    X_One_Tree_->GetDirectory()->cd();
    X_One_Tree_->Write();
	
    // 	estree_->GetDirectory()->cd();
    // 	estree_->Write();
    
    // 	mutree_->GetDirectory()->cd();
    // 	mutree_->Write();
	
    // 	//   etree_->GetDirectory()->cd();
    // 	//   etree_->Write();
	
    // 	trtree_->GetDirectory()->cd();
    // 	trtree_->Write();
	
    // 	xtree_->GetDirectory()->cd();
    // 	xtree_->Write();
	
    // 	jtree_->GetDirectory()->cd();
    // 	jtree_->Write();
	
    // 	if(doMC){
    // 		mctree_->GetDirectory()->cd();
    // 		mctree_->Write();
    // 	}
}//endjob



bool  JPsiPiPiPAT::isAbHadron(int pdgID) {
    
    if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
    return false;
    
}

bool 
JPsiPiPiPAT::isAMixedbHadron(int pdgID, int momPdgID) {
    
    if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) || 
        (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0)) 
        return true;
    return false;
    
}          

std::pair<int, float>   JPsiPiPiPAT::findCandMCInfo(reco::GenParticleRef genCand) {
    
    int momJpsiID = 0;
    float trueLife = -99.;
    //cout <<"externalmodule"<<endl;
    
    if (genCand->numberOfMothers()>0) {
        
        TVector3 trueVtx(0.0,0.0,0.0);
        TVector3 trueP(0.0,0.0,0.0);
        TVector3 trueVtxMom(0.0,0.0,0.0);
        
        trueVtx.SetXYZ(genCand->vertex().x(),genCand->vertex().y(),genCand->vertex().z());
        trueP.SetXYZ(genCand->momentum().x(),genCand->momentum().y(),genCand->momentum().z());
        
        bool aBhadron = false;
        reco::GenParticleRef Candmom = genCand->motherRef();       // find mothers
        if (Candmom.isNull()) {
            std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
            return result;
        } else {
            reco::GenParticleRef CandGrandMom = Candmom->motherRef();
            if (isAbHadron(Candmom->pdgId())) {
                if (CandGrandMom.isNonnull() && isAMixedbHadron(Candmom->pdgId(),CandGrandMom->pdgId())) {
                    momJpsiID = CandGrandMom->pdgId();
                    trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
                } else {
                    momJpsiID = Candmom->pdgId();
                    trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z());
                }
                aBhadron = true;
            } else {
                if (CandGrandMom.isNonnull() && isAbHadron(CandGrandMom->pdgId())) {
                    reco::GenParticleRef JpsiGrandgrandmom = CandGrandMom->motherRef();
                    if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(CandGrandMom->pdgId(),JpsiGrandgrandmom->pdgId())) {
                        momJpsiID = JpsiGrandgrandmom->pdgId();
                        trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
                    } else {
                        momJpsiID = CandGrandMom->pdgId();
                        trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
                    }
                    aBhadron = true;
                }
            }
            if (!aBhadron) {
                momJpsiID = Candmom->pdgId();
                trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z()); 
            }
        } 
        
        TVector3 vdiff = trueVtx - trueVtxMom;
        //trueLife = vdiff.Perp()*3.09688/trueP.Perp();
        trueLife = vdiff.Perp()*genCand->mass()/trueP.Perp();
    }
    std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
    return result;
    
}

double JPsiPiPiPAT::getSigmaOfLogdEdx(double logde)
{
	return 0.3;
}

float JPsiPiPiPAT::getEnergyLoss(const reco::TrackRef & track)
{
	if (iexception_dedx == 1)
		return 9999.;
	//   const DeDxDataValueMap & eloss = *energyLoss;
	const reco::DeDxDataValueMap & eloss = *energyLoss;
	return eloss[track].dEdx();
}

double JPsiPiPiPAT::nsigmaofdedx(const reco::TrackRef & track, double & theo,
		double & sigma)
{

	//   ch = (track->charge() > 0 ? 0 : 1);

	// no usable dE/dx if p > 2

	double nsigma = 99;
	if (iexception_dedx == 1)
		return nsigma;

	//   if(track->p() > 2) return nsigma;

	double m = 0.13957;
	double bg = track->p() / m;

	theo = getLogdEdx(bg);

	// !!!!!!
	int nhitr = track->numberOfValidHits();

	double meas = log(getEnergyLoss(track));
	sigma = getSigmaOfLogdEdx(theo) * pow(nhitr, -0.65);
	//   double errdedxTrk = eloss[trk1Ref].dEdxError();
	if (sigma > 0)
		nsigma = (meas - theo) / sigma;
	return nsigma;
}

double JPsiPiPiPAT::getLogdEdx(double bg)
{
	const double a = 3.25;
	const double b = 0.288;
	const double c = -0.852;

	double beta = bg / sqrt(bg * bg + 1);
	double dedx = log(a / (beta * beta) + b * log(bg) + c);

	return dedx;

}

double JPsiPiPiPAT::GetMass(const reco::TrackRef & track)
{
	double P = track->p();
	double C = 2.625;
	double K = 2.495;
	double I = getEnergyLoss(track);
	return sqrt((I - C) / K) * P;
}

double JPsiPiPiPAT::GetMassMC(const reco::TrackRef & track)
{
    double P = track->p();
    double C = 2.846;
    double K = 2.290;
    double I = getEnergyLoss(track);
    return sqrt((I - C) / K) * P;
}


//from UAMulti author:      Xavier Janssen

// void JPsiPiPiPAT::GetFwdGap(const edm::Event& iEvent , const edm::EventSetup& iSetup )
// {
//     using namespace std;
// 
// 
//     int nTowersHF_plus = 0;
//     int nTowersHF_minus = 0;
//     int nTowersHE_plus = 0;
//     int nTowersHE_minus = 0;
//     int nTowersHB_plus = 0;
//     int nTowersHB_minus = 0;
//     int nTowersEE_plus = 0;
//     int nTowersEE_minus = 0;
//     int nTowersEB_plus = 0;
//     int nTowersEB_minus = 0; 
//   //Sum(E)
//     double sumEHF_plus = 0.;
//     double sumEHF_minus = 0.;
//     double sumEHE_plus = 0.;
//     double sumEHE_minus = 0.;
//     double sumEHB_plus = 0.;
//     double sumEHB_minus = 0.;
//     double sumEEE_plus = 0.;
//     double sumEEE_minus = 0.;
//     double sumEEB_plus = 0.;
//     double sumEEB_minus = 0.;
//   // Sum(ET)
//     double sumETHF_plus = 0.;
//     double sumETHF_minus = 0.;
//     double sumETHE_plus = 0.;
//     double sumETHE_minus = 0.;
//     double sumETHB_plus = 0.;
//     double sumETHB_minus = 0.;
//     double sumETEE_plus = 0.;
//     double sumETEE_minus = 0.;
//     double sumETEB_plus = 0.;
//     double sumETEB_minus = 0.;
// 
//     
//     InputTag caloTowerTag_ = InputTag("towerMaker");
//     
//     const double energyThresholdHB_ = 1.5;
//     const double energyThresholdHE_ = 2.0;
//     const double energyThresholdHF_ = 4.0;
//     const double energyThresholdEB_ = 1.5;
//     const double energyThresholdEE_ = 2.5;
//     
//     
//   // Calo tower collection from event
//     edm::Handle<CaloTowerCollection> towerCollectionH;
//     iEvent.getByLabel(caloTowerTag_,towerCollectionH);
//     const CaloTowerCollection& towerCollection = *towerCollectionH;
// 
//   // Loop over calo towers
//     CaloTowerCollection::const_iterator calotower = towerCollection.begin();
//     CaloTowerCollection::const_iterator calotowers_end = towerCollection.end();
//     for(; calotower != calotowers_end; ++calotower) {
//         bool hasHCAL = false;
//         bool hasHF = false;
//         bool hasHE = false;
//         bool hasHB = false;
//         bool hasHO = false;
//         bool hasECAL = false;
//         bool hasEE = false;
//         bool hasEB = false;     
//         for(size_t iconst = 0; iconst < calotower->constituentsSize(); iconst++){
//             DetId detId = calotower->constituent(iconst);
//             if(detId.det()==DetId::Hcal){
//                 hasHCAL = true;
//                 HcalDetId hcalDetId(detId);
//                 if(hcalDetId.subdet()==HcalForward) hasHF = true;
//                 else if(hcalDetId.subdet()==HcalEndcap) hasHE = true;
//                 else if(hcalDetId.subdet()==HcalBarrel) hasHB = true;
//                 else if(hcalDetId.subdet()==HcalOuter) hasHO = true;  
//             } else if(detId.det()==DetId::Ecal){
//                 hasECAL = true;
//                 EcalSubdetector ecalSubDet = (EcalSubdetector)detId.subdetId();
//                 if(ecalSubDet == EcalEndcap) hasEE = true;
//                 else if(ecalSubDet == EcalBarrel) hasEB = true;
//             }
//         }
// 
//         int zside = calotower->zside();
//         double caloTowerEnergy = calotower->energy();
//      // FIXME
//      //double caloTowerET = calotower->et(primVtx.position());
//      //double caloTowerET = calotower->et(primVtx.z());
//         double caloTowerET = calotower->et();
// 
//      // HCAL: Towers made of at least one component from HB,HE,HF
//         if( hasHF && !hasHE ){
//             if( caloTowerEnergy >= energyThresholdHF_ ){
//                 if(zside >= 0){
//                     ++nTowersHF_plus;
//                     sumEHF_plus += caloTowerEnergy; 
//                     sumETHF_plus += caloTowerET;
//                 } else{
//                     ++nTowersHF_minus;
//                     sumEHF_minus += caloTowerEnergy;
//                     sumETHF_minus += caloTowerET;
//                 } 
//             }
//         } else if( hasHE && !hasHF && !hasHB ){
//             if( caloTowerEnergy >= energyThresholdHE_ ){
//                 if(zside >= 0){
//                     ++nTowersHE_plus;
//                     sumEHE_plus += caloTowerEnergy;
//                     sumETHE_plus += caloTowerET;
//                 } else{
//                     ++nTowersHE_minus;
//                     sumEHE_minus += caloTowerEnergy;
//                     sumETHE_minus += caloTowerET;
//                 }
//             }
//         } else if( hasHB && !hasHE ){
//             if( caloTowerEnergy >= energyThresholdHB_ ){
//                 if(zside >= 0){
//                     ++nTowersHB_plus;
//                     sumEHB_plus += caloTowerEnergy;
//                     sumETHB_plus += caloTowerET;
//                 } else{
//                     ++nTowersHB_minus;
//                     sumEHB_minus += caloTowerEnergy;
//                     sumETHB_minus += caloTowerET;
//                 }
//             }
//         }
// 
//      // ECAL: Towers made of at least one component from EB,EE
//         if( hasEE && !hasEB ){
//             if( caloTowerEnergy >= energyThresholdEE_ ){
//                 if(zside >= 0){
//                     ++nTowersEE_plus;
//                     sumEEE_plus += caloTowerEnergy;
//                     sumETEE_plus += caloTowerET;
//                 } else{
//                     ++nTowersEE_minus;
//                     sumEEE_minus += caloTowerEnergy;
//                     sumETEE_minus += caloTowerET;
//                 }
//             }
//         } else if( hasEB && !hasEE ){
//             if( caloTowerEnergy >= energyThresholdEB_ ){
//                 if(zside >= 0){
//                     ++nTowersEB_plus;
//                     sumEEB_plus += caloTowerEnergy;
//                     sumETEB_plus += caloTowerET;
//                 } else{
//                     ++nTowersEB_minus;
//                     sumEEB_minus += caloTowerEnergy;
//                     sumETEB_minus += caloTowerET;
//                 }
//             }
//         }
//     }
//     
//     CaloTowerInfo.Reset();
//     
//     CaloTowerInfo.eHfPos = sumEHF_plus;
//     CaloTowerInfo.eHfNeg = sumEHF_minus;
//     CaloTowerInfo.nHfTowersP = nTowersHF_plus;
//     CaloTowerInfo.nHfTowersN = nTowersHF_minus;
// }

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiPiPiPAT);

