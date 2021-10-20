#include <iostream>


// Dump data to the terminal
void dump( 	int identifier, int n,
						float* x, float* y, float* z,
						int* charge,
						float* px, float* py, float* pz ) {

	std::cout << identifier << std::endl
						<< *x					<< ' '
						<< *y					<< ' '
						<< *z					<< std::endl
						<< n					<< std::endl;

	// Loop over n particles
	for (int i=0; i<n; ++i) {
		std::cout << *(charge + i) 	<< ' '
							<< *(px + i) 			<< ' '
							<< *(py + i) 			<< ' '
							<< *(pz + i) 			<< std::endl;
	}


	return;
}
