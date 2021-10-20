#include <cmath>

#include "Event.h"

double mass( const Event& ev );


bool add( const Event& ev,
					float minMass, float maxMass,
					double& massSum, double& massSquare ) {

	double ptMass = mass(ev);

	if ( (ptMass >= minMass) && (ptMass <= maxMass) ) {
		massSum += ptMass;
		massSquare += std::pow( ptMass, 2 );
		return true;
	}

	else
		return false;
}
