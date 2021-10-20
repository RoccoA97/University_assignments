#include "Event.h"


Event::Event( int n, float x1, float y1, float z1 ):
	// initializations
 	eventId(n),
 	x(x1),
 	y(y1),
 	z(z1),
 	nPart(0) {
  	// allocate a buffer for particle pointers
  	pt = new Particle*[maxPart];
}


Event::~Event() {
  // delete all the particle pointers
  for ( unsigned int i = 0; i < nPart; ++i ) delete pt[i];

  // delete the pointers array
  delete[] pt;
}


void Event::add( float px, float py, float pz, int charge ) {

  // check for the number of particles, if maximum reached do nothing
	// and return
	if ( nPart >= maxPart ) return;

	// create the new particle and fill with data
	Particle* part = new Particle;
	part->px = px;
	part->py = py;
	part->pz = pz;
	part->charge = charge;

	// store the new particle pointer in the array and increase counter
	pt[nPart] = part;
	++nPart;

	return;
}


// get event id.
unsigned int Event::eventNumber() const {
	return eventId;
}


// get decay point coordinates
float Event::getX() const {
	return x;
}

float Event::getY() const {
	return y;
}

float Event::getZ() const {
	return z;
}


// get number of particles
unsigned int Event::nParticles() const {
	return nPart;
}


// get particle
const Event::Particle* Event::particle( unsigned int i ) const {

	if ( i < nPart )
    return pt[i];

	else
    return 0;
}
