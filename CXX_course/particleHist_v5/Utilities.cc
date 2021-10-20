#include "Utilities.h"
#include <cmath>

double Utilities::eFunc( double px, double py, double pz, double invMass ){
	return sqrt( pow(px,2) + pow(py,2) + pow(pz, 2) + pow(invMass,2) );
}

double Utilities::mFunc( double px, double py, double pz, double energy ) {

	double invMass2 = pow(energy,2) - ( pow(px,2) + pow(py,2) + pow(pz, 2) );

	if ( invMass2 > 0. )
		return sqrt( invMass2 );
	else
		return 0.;
}
