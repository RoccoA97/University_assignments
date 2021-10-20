#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

struct Event;

const Event* read( std::ifstream& iFile );
void dump( const Event& ev );
void clear( const Event* ev );
double eFunc( double px, double py, double pz, double invMass );
double mFunc( double px, double py, double pz, double energy );
double mass( const Event& ev );
bool add( const Event& ev, const float minMass, const float maxMass, double& massSum, double& massSquare );


int main( int argc, char* argv[] ) {

	const char* fileName = argv[1];
	std::ifstream iFile(fileName);

	unsigned int nEvent = 0;
	double massSum = 0;
	double massSquare = 0;
	double massMean;
	double massRMS;
	const float minMass = 0.490;
	const float maxMass = 0.505;

	const Event* ev;
	while ( (ev = read(iFile)) != 0 ) {
		// Comment the function dump(*ev) to prevent dump of data
		dump(*ev);
		if ( add( *ev, minMass, maxMass, massSum, massSquare ) ) ++nEvent;
		clear(ev);
	}

	massMean = massSum / nEvent;
	massRMS = std::sqrt( (massSquare/nEvent) - std::pow( (massSum/nEvent), 2. ) );

	std::cout << std::endl;
	std::cout << nEvent 	<< ' '
						<< massMean << ' '
						<< massRMS 	<< std::endl;


	return 0;
}
