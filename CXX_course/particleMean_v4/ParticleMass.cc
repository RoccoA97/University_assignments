#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "ParticleMass.h"
#include "Event.h"
#include "MassMean.h"

using namespace std;


ParticleMass::ParticleMass() {
}


ParticleMass::~ParticleMass() {
}


// function to be called at execution start
void ParticleMass::beginJob() {

  massVector.reserve( 2 );
  massVector.push_back( new MassMean( 0.490, 0.505 ) );
  massVector.push_back( new MassMean( 1.114, 1.118 ) );

  return;

}


// function to be called at execution end
void ParticleMass::endJob() {

  // loop over elements
  int n = massVector.size();
  int i;
  for ( i = 0; i < n; ++i ) {
    // get mass informations
    MassMean* mMean = massVector[i];
    // compute results
    mMean->compute();
    // get mean and rms masses and print results
    cout << i+1 << ' '  << mMean->getAcceptedEvent()  << ' '
                        << mMean->getMassMean()       << ' '
                        << mMean->getMassRMS()        << endl;
  }

  return;

}


// function to be called for each event
void ParticleMass::process( const Event& ev ) {

  // loop over masses and pass event to each of them
  unsigned int n = massVector.size();
  unsigned int i;
  for ( i = 0; i < n; ++i ) massVector[i]->add( ev );

  return;

}
