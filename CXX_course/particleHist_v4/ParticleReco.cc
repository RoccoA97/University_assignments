#include "ParticleReco.h"
#include "Event.h"
#include "Utilities.h"
#include "Constants.h"
#include <iostream>
#include <math.h>

using namespace std;


ParticleReco::ParticleReco() {
}


ParticleReco::~ParticleReco() {
}


// recompute tag informations for new event
void ParticleReco::update( const Event& ev ) {

  // set default quantities
  type   = unknown;
  energy = -1.0;
  mass   = -1.0;

  // code to compute total energy and invariant mass for the two
  // mass hypotheses for the decay products

  // compare invariant masses with known values and set final values
  // ( type, energy and mass )

  // variables to loop over particles
	// positive / negative track counters
  int nCharge = 0;
	int pCharge = 0;

  // variables for momentum sums
	double pxSum = 0.;
	double pySum = 0.;
	double pzSum = 0.;

  // variables for energy sums, for K0 and Lambda0
	double K0Energy = 0.;
	double L0Energy = 0.;
	double K0Mass = 0.;
	double L0Mass = 0.;

  // loop over particles - momenta
	for ( unsigned int i = 0; i < ev.nParticles(); ++i ) {
    // get particle pointer
		const Event::Particle* part = ev.particle(i);

		// update momentum sums
		pxSum += (part->px);
		pySum += (part->py);
		pzSum += (part->pz);

    // update energy sums, for K0 and Lambda0:
    // pion mass for negative particle,
  	// pion or proton mass for positive particle,
    // for K0 or Lambda0 respectively
		K0Energy += Utilities::eFunc( (part->px), (part->py), (part->pz), Constants::massPion );

    if ( (part->charge) < 0 )
			L0Energy += Utilities::eFunc( (part->px), (part->py), (part->pz), Constants::massPion );
		if ( (part->charge) > 0 )
			L0Energy += Utilities::eFunc( (part->px), (part->py), (part->pz), Constants::massProton );


		// update positive/negative track counters
		if ( (part->charge) < 0 ) ++nCharge;
		if ( (part->charge) > 0 ) ++pCharge;
	}

  // check for exactly one positive and one negative track
	// otherwise return negative (unphysical) invariant mass
	if ( nCharge==1 && pCharge==1 ) {
		K0Mass = Utilities::mFunc( pxSum, pySum, pzSum, K0Energy );
		L0Mass = Utilities::mFunc( pxSum, pySum, pzSum, L0Energy );

		double K0Diff = fabs( Constants::massK0 - K0Mass );
		double L0Diff = fabs( Constants::massLambda0 - L0Mass );

		if ( K0Diff < L0Diff ) {
      type = K0;
      energy = K0Energy;
      mass = K0Mass;
    }
		else {
      type = Lambda0;
      energy = L0Energy;
      mass = L0Mass;
    }
	}

  return;

}


float ParticleReco::totalEnergy() {
  // check for new event and return result
  check();
  return energy;
}


float ParticleReco::invariantMass() {
  // check for new event and return result
  check();
  return mass;
}


ParticleReco::ParticleType ParticleReco::particleType() {
  // check for new event and return result
  check();
  return type;
}
