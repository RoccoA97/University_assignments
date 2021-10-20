#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "ParticleMass.h"
#include "Event.h"
#include "MassMean.h"

#include "TH1F.h"
#include "TFile.h"

using namespace std;

double mass( const Event& ev );


ParticleMass::ParticleMass() {
}


ParticleMass::~ParticleMass() {
}


// function to be called at execution start
void ParticleMass::beginJob() {

  massVector.reserve( 2 );
  pCreate( "K0Hist", 0.495, 0.500 );
  pCreate( "L0Hist", 1.115, 1.116 );

  return;

}


// function to be called at execution end
void ParticleMass::endJob() {

  // save current working area
  TDirectory* currentDir = gDirectory;
  // open histogram file
  TFile* file = new TFile( "hist.root", "CREATE" );


  // loop over elements
  int n = massVector.size();
  int i;
  for ( i = 0; i < n; ++i ) {
    // get mass informations
    Particle* mMean = massVector[i];
    // compute results
    mMean->massMeanPtr->compute();
    // get mean and rms masses and print results
    cout << i+1 << ' '  << mMean->massMeanPtr->getAcceptedEvent()  << ' '
                        << mMean->massMeanPtr->getMassMean()       << ' '
                        << mMean->massMeanPtr->getMassRMS()        << endl;

    mMean->rootHist->Write();
  }

  file->Close();
  delete file;
  currentDir->cd();

  return;

}


// function to be called for each event
void ParticleMass::process( const Event& ev ) {
  
  // loop over masses and pass event to each of them
  unsigned int n = massVector.size();
  unsigned int i;
  for ( i = 0; i < n; ++i ) {
    if ( massVector[i]->massMeanPtr->add( ev ) )
      massVector[i]->rootHist->Fill( mass( ev ) );
  }

  return;

}


void ParticleMass::pCreate( const std::string& name, float min, float max ) {

  const char* hName = name.c_str();

  unsigned int nBin = 100;

  Particle* part = new Particle;
  part->name = name;
  part->massMeanPtr = new MassMean( min, max );
  part->rootHist = new TH1F( hName, hName, nBin, min, max );
  massVector.push_back( part );

  return;

}
