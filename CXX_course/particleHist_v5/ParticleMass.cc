#include "ParticleMass.h"

#include "Event.h"
#include "AnalysisInfo.h"
#include "AnalysisFactory.h"
#include "MassMean.h"
#include "ParticleReco.h"

#include "util/include/ActiveObserver.h"

#include "TH1F.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

double mass( const Event& ev );


class ParticleMassFactory: public AnalysisFactory::AbsFactory {
 public:
  // assign "hist" as name for this analyzer and factory
  ParticleMassFactory(): AnalysisFactory::AbsFactory( "hist" ) {}
  // create an ParticleMass when builder is run
  virtual AnalysisSteering* create( const AnalysisInfo* info ) {
    return new ParticleMass( info );
  }
};
// create a global ParticleMassFactory, so that it is created and registered
// before main execution start:
// when the AnalysisFactory::create function is run,
// a ParticleMassFactory will be available with name "hist".
static ParticleMassFactory er;


ParticleMass::ParticleMass( const AnalysisInfo* info ):
  AnalysisSteering( info ) {
}


ParticleMass::~ParticleMass() {
}


// function to be called at execution start
void ParticleMass::beginJob() {

  // create mass distributions for different mass ranges
  ifstream file( aInfo->value( "ranges" ).c_str() );
  string name;
  float minMass;
  float maxMass;

  massVector.reserve( 2 );
  while ( file >> name >> minMass >> maxMass ) pCreate( name, minMass, maxMass );

  return;

}


// function to be called at execution end
void ParticleMass::endJob() {

  // save current working area
  TDirectory* currentDir = gDirectory;
  // open histogram file
  TFile* file = new TFile( aInfo->value( "hist" ).c_str(), "CREATE" );


  // loop over elements
  int n = massVector.size();
  int i;
  for ( i = 0; i < n; ++i ) {
    // get mass informations
    Particle* mMean = massVector[i];
    // compute results
    mMean->massMeanPtr->compute();
    // get mean and rms masses and print results
    cout << i << " " << mMean->massMeanPtr->getMassMean() << " " << mMean->massMeanPtr->getMassRMS() << endl;
    // print number of events
    cout << mMean->massMeanPtr->getAcceptedEvent() << endl;

    mMean->rootHist->Write();
  }

  file->Close();
  delete file;
  currentDir->cd();

  return;

}


// function to be called for each event
void ParticleMass::update( const Event& ev ) {
  // loop over masses and pass event to each of them

  static ParticleReco* pRec = ParticleReco::instance();
  float partMass = pRec->invariantMass();

  unsigned int n = massVector.size();
  unsigned int i;
  for ( i = 0; i < n; ++i ) {
    massVector[i]->massMeanPtr->add( ev );

    if ( partMass > 0)
    massVector[i]->rootHist->Fill( partMass );
  }

  return;
}

void ParticleMass::pCreate( const std::string& name, float min, float max ) {

  const char* hName = ("mass" + name).c_str();

  Particle* part = new Particle;
  part->name = name;
  part->massMeanPtr = new MassMean( min, max );
  part->rootHist = new TH1F( hName, hName, 100, min, max );
  massVector.push_back( part );

  return;

}
