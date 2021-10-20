#include <cmath>

#include "Event.h"


// compute energy from momentum x,y,z components and invariant mass
double eFunc( double px, double py, double pz, double invMass ) {

	return sqrt( pow(px,2) + pow(py,2) + pow(pz, 2) + pow(invMass,2) );
}


// compute invariant mass from momentum x,y,z components and energy
double mFunc( double px, double py, double pz, double energy ) {

	double invMass2 = std::pow(energy,2) - ( std::pow(px,2) + std::pow(py,2) + std::pow(pz, 2) );

	if ( invMass2 > 0. )
		return sqrt( invMass2 );

	else
		return 0.;
}


const double massPion    = 0.1395706;   // GeV/c^2
const double massProton  = 0.938272;    // GeV/c^2
const double massK0      = 0.497611;    // GeV/c^2
const double massLambda0 = 1.115683;    // GeV/c^2


double mass( const Event& ev ) {

  // positive / negative track counters
	int nCharge = 0;
	int pCharge = 0;

  // variables for momentum sums
	double pxSum = 0;
	double pySum = 0;
	double pzSum = 0;

  // variables for energy sums, for K0 and Lambda0
	double K0Energy;
	double L0Energy;
	double K0Mass;
	double L0Mass;

  // loop over particles - momenta
	for ( int i = 0; i < ev.nPart; ++i ) {
    // get particle pointer
		const Particle* part = ev.pt[i];

    // update momentum sums
		pxSum += (part->px);
		pySum += (part->py);
		pzSum += (part->pz);

    // update energy sums, for K0 and Lambda0:
		// pion mass for negative particle,
    // pion or proton mass for positive particle,
    // for K0 or Lambda0 respectively
		K0Energy += eFunc( (part->px), (part->py), (part->pz), massPion );

		if ( (part->charge) < 0 )
			L0Energy += eFunc( (part->px), (part->py), (part->pz), massPion );
		if ( (part->charge) > 0 )
			L0Energy += eFunc( (part->px), (part->py), (part->pz), massProton );

    // update positive/negative track counters
		if ( (part->charge) < 0 ) ++nCharge;
		if ( (part->charge) > 0 ) ++pCharge;
	}

  // check for exactly one positive and one negative track
  // otherwise return negative (unphysical) invariant mass
	if ( nCharge==1 && pCharge==1 ) {
		K0Mass = mFunc( pxSum, pySum, pzSum, K0Energy );
		L0Mass = mFunc( pxSum, pySum, pzSum, L0Energy );

		float K0Diff = fabs( massK0 - K0Mass );
		float L0Diff = fabs( massLambda0 - L0Mass );

		if ( K0Diff < L0Diff )
			return K0Mass;

		else
			return L0Mass;
	}


	else return -1;
}
