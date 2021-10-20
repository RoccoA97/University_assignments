#include <fstream>


// Read data from file
int read( std::ifstream& iFile,
					float* x, float* y, float* z,
					int* charge,
					float* px, float* py, float* pz ) {

	int n;

	iFile >> *x >> *y >> *z >> n;

	for (int i = 0; i < n; ++i) {
		iFile >> *(charge + i) >> *(px + i) >> *(py + i) >> *(pz + i);
	}

	return n;
}
