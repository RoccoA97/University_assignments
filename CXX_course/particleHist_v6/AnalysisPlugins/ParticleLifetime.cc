#include "ParticleLifetime.h"

#include "../AnalysisFramework/Event.h"
#include "../AnalysisFramework/AnalysisInfo.h"
#include "../AnalysisFramework/AnalysisFactory.h"
#include "../AnalysisObjects/MassMean.h"
#include "../AnalysisObjects/ParticleReco.h"
#include "../AnalysisObjects/LifetimeFit.h"
#include "../util/include/ActiveObserver.h"

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
  // assign "hist" as name for this analyzer and factory
  ParticleLifetimeFactory(): AnalysisFactory::AbsFactory( "time" ) {}
  // create an ParticleLifetime when builder is run
  virtual AnalysisSteering* create( const AnalysisInfo* info ) {
    return new ParticleLifetime( info );
  }
};
// create a global ElementRecoFactory, so that it is created and registered
// before main execution start:
// when the AnalysisFactory::create function is run,
// an ElementRecoFactory will be available with name "plot".
static ParticleLifetimeFactory er;


ParticleLifetime::ParticleLifetime( const AnalysisInfo* info ):
  AnalysisSteering( info ) {
}


ParticleLifetime::~ParticleLifetime() {
}


// function to be called at execution start
void ParticleLifetime::beginJob() {

  // create mass distributions for different mass ranges
  ifstream file( "particleFitters" );
  string name;
  double minMass;
  double maxMass;
  double minTime;
  double maxTime;
  double minScan;
  double maxScan;
  double scanStep;

  timeVector.reserve( 2 );
  while ( file >> name >> minMass >> maxMass >> minTime >> maxTime >> minScan >> maxScan >> scanStep )
    pCreate( name, minMass, maxMass, minTime, maxTime, minScan, maxScan, scanStep );

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
    // print number of events, mean life time and error
    cout << mMean->lifetimeFitPtr->nEvents()        << ' '
         << mMean->lifetimeFitPtr->lifeTime()       << ' '
         << mMean->lifetimeFitPtr->lifeTimeError()  << endl;

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

void ParticleLifetime::pCreate( const std::string& name,
                                double minMass, double maxMass,
                                double minTime, double maxTime,
                                double minScan, double maxScan, double scanStep ) {

  const char* hName = ("time" + name).c_str();

  Particle* part = new Particle;
  part->name = name;
  part->lifetimeFitPtr = new LifetimeFit( minMass, maxMass, minTime, maxTime, minScan, maxScan, scanStep );
  part->rootHist = new TH1F( hName, hName, 100, minTime, maxTime );
  timeVector.push_back( part );

  return;

}
