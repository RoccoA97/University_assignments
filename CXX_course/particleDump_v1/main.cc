#include <fstream>


int read( std::ifstream& file,
					float* x, float* y, float* z,
					int* charge,
					float* px, float* py, float* pz );

void dump( 	int identifier, int n,
						float* x, float* y, float* z,
						int* charge,
						float* px, float* py, float* pz );


int main( int argc, char* argv[] ) {

	const char* fileName = argv[1];
	std::ifstream iFile(fileName);

	// Static arrays and pointers to contain data from file
	// x,y,z 				decay point coordinates
	// px,py,pz 		momenta components
	// charge				electric charges of particles
	// eventNumber	event Id
	// nPart				number of particles in a single event
	float* x = new float;
	float* y = new float;
	float* z = new float;
	float px[10];
	float py[10];
	float pz[10];
	int charge[10];
	int eventNumber;
	int nPart;

	// Read from file and dump
	while ( iFile >> eventNumber  ) {
		nPart = read( iFile, x, y, z, charge, px, py, pz );
		dump( eventNumber, nPart, x, y, z, charge, px, py, pz );
	}


	return 0;
}
