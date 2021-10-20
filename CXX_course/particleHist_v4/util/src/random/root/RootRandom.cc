#include "RootRandom.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

static RootRandom rr;

RootRandom::RootRandom() {
  if ( verbosity >= 1 ) cout << "create RootRandom " << this << endl; 
  rand = new TRandom;
}


RootRandom::~RootRandom() {
}


// implementation of random number generation functions

void RootRandom::set( unsigned int seed ) {
  if ( verbosity >= 2 ) cout << "RootRandom::set called" << endl;
  rand->SetSeed( seed );
  return;
}


float RootRandom::generate( Random::probability p, float a, float b ) {
  if ( verbosity >= 2 ) cout << "RootRandom::generate called" << endl;
  if ( rand == 0 ) rand = new TRandom;
  switch ( p ) {
  case Random::Flat:
    return rand->Uniform( a, b );
  case Random::Gauss:
    return rand->Gaus( a, b );
  }
  return 0.0;
}


