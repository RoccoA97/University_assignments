#include "Random.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

int Random::verbosity = 0;

Random::Random() {
  if ( verbosity >= 1 ) cout << "create Random " << this << endl; 
  Random*& i = instRef();
  if ( i == 0 ) i = this;
}


Random::~Random() {
}


// random number generation functions

void Random::setSeed( unsigned int seed ) {
  if ( verbosity >= 3 ) cout << "Random::setSeed " << endl; 
  // get concrete random generator pointer
  // save a static reference to avoid repeating function call
  static Random* i = instance();
  if ( verbosity >= 3 ) cout << "Random::setSeed " << i << endl; 
  if ( i == 0 ) {
    cout << "concrete random generator not set" << endl;
  }
  i->set( seed );
  return;
}


float Random::flat( double min , double max ) {
  if ( verbosity >= 3 ) cout << "Random::flat " << endl; 
  // get concrete random generator pointer
  // save a static reference to avoid repeating function call
  static Random* i = instance();
  if ( verbosity >= 3 ) cout << "Random::flat " << i << endl; 
  if ( i == 0 ) {
    cout << "concrete random generator not set" << endl;
    return ( min + max ) / 2.0;
  }
  return i->generate( Flat, min, max );
}


float Random::gauss( double mean, double rms ) {
  if ( verbosity >= 3 ) cout << "Random::gauss " << endl; 
  // get concrete random generator pointer
  // save a static reference to avoid repeating function call
  static Random* i = instance();
  if ( verbosity >= 3 ) cout << "Random::gauss " << i << endl; 
  if ( i == 0 ) {
    cout << "concrete random generator not set" << endl;
    return mean;
  }
  return i->generate( Gauss, mean, rms );
}


// implementation of random number generation functions

void Random::set( unsigned int seed ) {
  if ( verbosity >= 2 ) cout << "StdRandom::set called" << endl;
  srandom( seed );
  return;
}


float Random::generate( Random::probability p, float a, float b ) {
  if ( verbosity >= 2 ) cout << "StdRandom::generate called" << endl;
  float x = random() * 1.0 / RAND_MAX;
  switch ( p ) {
  case Random::Flat:
    return a + ( x * ( b - a ) );
  case Random::Gauss:
    float y = random() * 1.0 / RAND_MAX;
    return a - ( b * sqrt( -log( 1.0 - x ) ) * cos( 2.0 * M_PI * y ) );
  }
  return 0.0;
}


// get a pointer to a random generator,
// eventually to be assigned a concrete random generator
Random* Random::instance() {
  static Random*& obj = instRef();
  if ( obj == 0 ) new Random;
  return obj;
}


// get a pointer to a random generator,
// eventually to be assigned a concrete random generator
Random*& Random::instRef() {
  if ( verbosity >= 1 ) cout << "Random::instance " << endl; 
  // the pointer is created only once, the first time "instance()" is called
  static Random* p = 0;
  return p;
}

