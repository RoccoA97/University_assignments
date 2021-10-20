#ifndef RootRandom_H
#define RootRandom_H

#include "Random.h"
#include "TRandom.h"

class RootRandom: public Random {

 public:

  RootRandom();
  ~RootRandom();

 protected:

  // redeclaration of random number generation functions
  virtual void set( unsigned int seed );
  virtual float generate( Random::probability p, float a, float b );

 private:

  // actual ROOT generator
  TRandom* rand;

  RootRandom           ( const RootRandom& x );
  RootRandom& operator=( const RootRandom& x );

};

#endif

