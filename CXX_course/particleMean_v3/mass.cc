#include <cmath>

#include "Event.h"
#include "MassMean.h"
#include "Utilities.h"
#include "Constants.h"


double mass( const Event& ev ) {

  // variables to loop over particles
  // positive / negative track counters
	int nCharge = 0;
	int pCharge = 0;

  // variables for momentum sums
	double pxSum = 0;
	double pySum = 0;
	double pzSum = 0;

  // variables for energy sums, for K0 and Lambda0
	double K0Energy = 0;
	double L0Energy = 0;
	double K0Mass = 0;
	double L0Mass = 0;

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
		double K0Diff = std::fabs( Constants::massK0 - K0Mass );
		double L0Diff = std::fabs( Constants::massLambda0 - L0Mass );

		if ( K0Diff < L0Diff )
			return K0Mass;

		else
			return L0Mass;
		}

	else
		return -1;
}
