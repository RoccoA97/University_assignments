#include <iostream>
#include <fstream>

struct Event;

const Event* read( std::ifstream& iFile );
void dump( const Event& ev );
void clear( const Event* ev );


int main( int argc, char* argv[] ) {

	const char* fileName = argv[1];
	std::ifstream iFile(fileName);

	// Read from file and dump
	const Event* ev;
	while ( (ev = read(iFile)) != 0 ) {
		dump(*ev);
		clear(ev);
	}


	return 0;
}
