#include <iostream>
#include <fstream>

#include "Event.h"
#include "MassMean.h"
#include "Utilities.h"
#include "Constants.h"

const Event* read( std::ifstream& iFile );
void dump( const Event& ev );
double mass( const Event& ev );


int main( int argc, char* argv[] ) {

	const char* fileName = argv[1];
	std::ifstream iFile(fileName);

	MassMean K0(0.490, 0.505);
	MassMean L0(1.114, 1.118);

	const Event* ev;
	while ( (ev = read(iFile)) != 0 ) {
		// Comment the function dump(*ev) to prevent dump of data
		K0.add(*ev);
		L0.add(*ev);
		dump(*ev);
		delete ev;
	}

	K0.compute();
	L0.compute();

	std::cout << std::endl;
	std::cout << K0.getAcceptedEvent() << ' ' << K0.getMassMean() << ' ' << K0.getMassRMS() << std::endl;
	std::cout << L0.getAcceptedEvent() << ' ' << L0.getMassMean() << ' ' << L0.getMassRMS() << std::endl;

	return 0;
}
