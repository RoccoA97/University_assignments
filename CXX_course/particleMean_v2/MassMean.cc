#include "MassMean.h"
#include "Event.h"

#include <cmath>

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

  // return true or false for invariant mass inside or outside the range
  double partMass = mass(ev);
  if ( (partMass >= minMass) && (partMass <= maxMass) ) {

  	massSum += partMass;
  	massSquare += std::pow( partMass, 2. );
  	++acceptedEvent;

  	return true;
  }

  else
		return false;
}


// compute mean and rms
void MassMean::compute() {
	massMean = massSum / acceptedEvent;
	massRMS = std::sqrt( (massSquare/acceptedEvent) - std::pow( (massSum/acceptedEvent), 2. ) );
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
