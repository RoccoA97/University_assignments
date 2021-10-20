#include <cmath>

#include "Utilities.h"


double Utilities::eFunc( double px, double py, double pz, double invMass ){
	return std::sqrt( std::pow(px,2) + std::pow(py,2) + std::pow(pz, 2) + std::pow(invMass,2) );
}


double Utilities::mFunc( double px, double py, double pz, double energy ) {

	double invMass2 = std::pow(energy,2) - ( std::pow(px,2) + std::pow(py,2) + std::pow(pz, 2) );

	if ( invMass2 > 0. )
		return std::sqrt( invMass2 );

	else
		return 0.;
}
