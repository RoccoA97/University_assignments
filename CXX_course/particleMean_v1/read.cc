#include <fstream>

#include "Event.h"


// Read from file
const Event* read( std::ifstream& iFile ) {

	// Variable to get values
	float f;
	Event* ev;

	if ( (iFile >> f) != 0 )
		ev = new Event;

	else
		return 0;

	// Fill ev members
	ev->eventId = static_cast<int>(f);
	iFile >> f;
	ev->x = f;
	iFile >> f;
	ev->y = f;
	iFile >> f;
	ev->z = f;
	iFile >> f;
	ev->nPart = static_cast<int>(f);

	// Create a dynamic array of pointers to Particle
	ev->pt = new Particle*[ev->nPart];

	for ( int j = 0; j < ev->nPart; ++j ) {
		// Create a pointer to Particle
		// then assign it to the object "part"
		// to get data from file
		Particle* part = ev->pt[j] = new Particle;
		iFile >> f;
		part->charge = static_cast<int>(f);
		iFile >> f;
		part->px = f;
		iFile >> f;
		part->py = f;
		iFile >> f;
		part->pz = f;
	}


	return ev;
}
