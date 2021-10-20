#include "ParticleLifetime.h"

#include "Event.h"
#include "AnalysisInfo.h"
#include "AnalysisFactory.h"
#include "MassMean.h"
#include "ParticleReco.h"
#include "LifetimeFit.h"

#include "util/include/ActiveObserver.h"

#include "TH1F.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

double mass( const Event& ev );


class ParticleLifetimeFactory: public AnalysisFactory::AbsFactory {
 public:
  // assign "time" as name for this analyzer and factory
  ParticleLifetimeFactory(): AnalysisFactory::AbsFactory( "time" ) {}
  // create an ParticleLifetime when builder is run
  virtual AnalysisSteering* create( const AnalysisInfo* info ) {
    return new ParticleLifetime( info );
  }
};
// create a global ParticleLifetimeFactory, so that it is created and registered
// before main execution start:
// when the AnalysisFactory::create function is run,
// a ParticleLifetimeFactory will be available with name "time".
static ParticleLifetimeFactory er;


ParticleLifetime::ParticleLifetime( const AnalysisInfo* info ):
  AnalysisSteering( info ) {
}


ParticleLifetime::~ParticleLifetime() {
}


// function to be called at execution start
void ParticleLifetime::beginJob() {

  timeVector.reserve( 2 );
  pCreate( "K0Hist", 0.490, 0.505, 10.0, 500.0  );
  pCreate( "L0Hist", 1.114, 1.118, 10.0, 1000.0 );

  return;

}


// function to be called at execution end
void ParticleLifetime::endJob() {

  // save current working area
  TDirectory* currentDir = gDirectory;
  // open histogram file
  TFile* file = new TFile( aInfo->value( "time" ).c_str(), "CREATE" );


  // loop over elements
  int n = timeVector.size();
  int i;
  for ( i = 0; i < n; ++i ) {
    // get mass informations
    Particle* mMean = timeVector[i];
    // compute results
    mMean->lifetimeFitPtr->compute();
    // print number of events
    cout << mMean->lifetimeFitPtr->getAcceptedEvent() << endl;

    mMean->rootHist->Write();
  }



  file->Close();
  delete file;
  currentDir->cd();

  return;

}


// function to be called for each event
void ParticleLifetime::update( const Event& ev ) {

  // loop over masses and pass event to each of them

  static ProperTime* pTime = ProperTime::instance();
  float partTime = pTime->decayTime();

  unsigned int n = timeVector.size();
  unsigned int i;
  for ( i = 0; i < n; ++i ) {
    // control condition for masses
    if ( timeVector[i]->lifetimeFitPtr->add( ev ) )
      if ( partTime > 0)
        timeVector[i]->rootHist->Fill( partTime );
  }

  return;

}


void ParticleLifetime::pCreate( const std::string& name, float min, float max, float timeMin, float timeMax ) {

  const char* hName = ("time" + name).c_str();

  Particle* part = new Particle;
  part->name = name;
  part->lifetimeFitPtr = new LifetimeFit( min, max );
  part->rootHist = new TH1F( hName, hName, 100, timeMin, timeMax );
  timeVector.push_back( part );

  return;

}
