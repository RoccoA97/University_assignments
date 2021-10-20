#include "Event.h"

Event::Event( int n, float x1, float y1, float z1 ):
  // initializations
  eventId( n ),
  x ( x1 ),
  y ( y1 ),
  z ( z1 ) {
  // allocate a buffer for particle pointers
  pt.reserve( 10 );
}


Event::~Event() {
  unsigned int i;
  // delete all the particle pointers
  for ( i = 0; i < pt.size(); ++i ) delete pt[i];
}


void Event::add( float px, float py, float pz, int charge ) {

  // create the new particle and fill with data
  Particle* part = new Particle;
  part->charge = charge;
  part->px     = px;
  part->py     = py;
  part->pz     = pz;

  // store the new particle pointer in the array and increase counter
  pt.push_back( part );

  return;

}


// get event id.
int Event::eventNumber() const {
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
  return pt.size();
}


// get particle
const Event::Particle* Event::particle( unsigned int i ) const {
  // check for required particle being inside the array:
  // - if required particle inside the array return its pointer
  // - if not return null
  if ( i < pt.size() ) return pt[i];
  return 0;
}
