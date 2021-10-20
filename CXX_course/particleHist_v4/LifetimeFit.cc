#include "ProperTime.h"
#include "LifetimeFit.h"
#include "Event.h"
#include "ParticleReco.h"

#include <cmath>
#include <iostream> //

double mass( const Event& ev );


LifetimeFit::LifetimeFit( float min, float max ):
	minMass(min),
	maxMass(max),
	acceptedEvent(0) {}


LifetimeFit::~LifetimeFit() {}


// add data from a new event
bool LifetimeFit::add( const Event& ev ) {

	static ParticleReco* pRec = ParticleReco::instance();
  float partMass = pRec->invariantMass();

  if ( (partMass >= minMass) && (partMass <= maxMass) ) {
  	++acceptedEvent;
  	return true;
  }

  else return false;

}


// compute
void LifetimeFit::compute() {}


// return number of selected events
unsigned int LifetimeFit::getAcceptedEvent() const {
	return acceptedEvent;
}
