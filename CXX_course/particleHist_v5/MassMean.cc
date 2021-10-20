#include "MassMean.h"
#include "Event.h"
#include "ParticleReco.h"

#include <cmath>
#include <iostream> //

double mass( const Event& ev );


MassMean::MassMean( float min, float max ):
	minMass(min),
	maxMass(max),
	acceptedEvent(0),
	massSum(0),
	massSquare(0) {}


MassMean::~MassMean() {}


// add data from a new event
bool MassMean::add( const Event& ev ) {

	static ParticleReco* pRec = ParticleReco::instance();
  float partMass = pRec->invariantMass();

  if ( (partMass >= minMass) && (partMass <= maxMass) ) {

  	massSum += partMass;
  	massSquare += pow( partMass, 2 );
  	++acceptedEvent;

  	return true;
  }

  else return false;

}


// compute mean and rms
void MassMean::compute() {
	massMean = massSum / acceptedEvent;
	massRMS = sqrt( (massSquare/acceptedEvent) - pow( (massSum/acceptedEvent), 2 ) );
}


// return number of selected events
unsigned int MassMean::getAcceptedEvent() const {
	return acceptedEvent;
}


// return mean and rms
double MassMean::getMassMean() const {
	return massMean;
}

double MassMean::getMassRMS() const {
	return massRMS;
}
