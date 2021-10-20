#include "ProperTime.h"
#include "LifetimeFit.h"
#include "ParticleReco.h"

#include "../AnalysisFramework/Event.h"
#include "../AnalysisUtilities/QuadraticFitter.h"

#include <cmath>

double mass( const Event& ev );


LifetimeFit::LifetimeFit( float min, float max,
 													double minTime, double maxTime,
													double minScan, double maxScan, double scanStep ):
													minMass(min),
													maxMass(max),
													minTime(minTime),
													maxTime(maxTime),
													minScan(minScan),
													maxScan(maxScan),
													scanStep(scanStep) {}


LifetimeFit::~LifetimeFit() {}


// add data from a new event
bool LifetimeFit::add( const Event& ev ) {

	static ProperTime* ptProperTime = ProperTime::instance();
	static ParticleReco* pRec = ParticleReco::instance();
  float partMass = pRec->invariantMass();
	double partDecayTime = ptProperTime->decayTime();

  if ( (partMass >= minMass) && (partMass <= maxMass) ) {
    if ( (partDecayTime >= minTime) && (partDecayTime <= maxTime) ) {
      decayTimes.push_back( partDecayTime );
      return true;
    }

    else
      return false;
  }

  else
    return false;

}


// compute
void LifetimeFit::compute() {

	QuadraticFitter qFit;

	double t = minScan;

	// add data to fit
	while ( t <= maxScan ) {
    double L = 0;
		for ( unsigned int i=0; i < decayTimes.size(); ++i ) {
			L += ( ( decayTimes[i] / t ) + log( t ) + log( exp( - minTime / t ) - exp( - maxTime / t )));
		}

		qFit.add(t,L);
		t += scanStep;
	}

	ptLifeTimeMean = - qFit.b() / ( 2. * qFit.c() );
	ptLifeTimeError = 1. / sqrt( 2. * qFit.c() );

	qFit.clear();
	qFit.reset();

	return;

}


// return mean life time
double LifetimeFit::lifeTime() const {
	return ptLifeTimeMean;
}


// return life time error
double LifetimeFit::lifeTimeError() const {
	return ptLifeTimeError;
}


// return number of selected events
unsigned int LifetimeFit::nEvents() const {
	return decayTimes.size();
}
