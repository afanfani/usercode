// -*- C++ -*-
//
// Package:    FakeTestAnalyzer
// Class:      FakeTestAnalyzer
// 
/**\class FakeTestAnalyzer FakeTestAnalyzer.cc AFanfani/FakeTestAnalyzer/src/FakeTestAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  
//         Created:  Mon Oct  6 11:03:20 CEST 2008
// $Id$
//
//


// system include files
#include <memory>

#include <cstring>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace std;

//
// class decleration
//

class FakeTestAnalyzer : public edm::EDAnalyzer {
   public:
      explicit FakeTestAnalyzer(const edm::ParameterSet&);
      ~FakeTestAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
     int sleep_time;            // integer for the time to sleep
     double output_size;        // for text output, decide the size (0.1 MB a unit)
     std::string text_filename; // the filename for the output file

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
FakeTestAnalyzer::FakeTestAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  text_filename = iConfig.getUntrackedParameter<string>("OutTextFileName");
  output_size   = iConfig.getUntrackedParameter<double>("OutputSize");
  sleep_time    = iConfig.getUntrackedParameter<int>("SleepTime");

  cout<<"OutTextFileName      "<<text_filename<<endl;
  cout<<"OutputSize        "<<output_size<<endl;
  cout<<"SleepTime         "<<sleep_time<<endl;

}


FakeTestAnalyzer::~FakeTestAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  int i, j;
  char * filename = new char[text_filename.length() + 1];
  strcpy(filename, text_filename.c_str());

  ofstream myfile;
  myfile.open (filename);
  for (i=0; i< output_size/0.1; i++){
      for (j=0; j<3303; j++){
        myfile << "This is a grid test file, " <<output_size<<"MB"<<endl;
      }
  }
  myfile.close();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
FakeTestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // sleep time
   sleep(sleep_time);

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
FakeTestAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FakeTestAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(FakeTestAnalyzer);
