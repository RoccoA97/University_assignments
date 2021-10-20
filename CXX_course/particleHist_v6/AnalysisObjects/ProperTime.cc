#include "ProperTime.h"
#include "ParticleReco.h"

#include "../AnalysisFramework/Event.h"
#include "../AnalysisUtilities/Utilities.h"
#include "../AnalysisUtilities/Constants.h"

#include <iostream>
#include <math.h>

using namespace std;


ProperTime::ProperTime() {
}


ProperTime::~ProperTime() {
}


// recompute tag informations for new event
void ProperTime::update( const Event& ev ) {

  // set default quantities
  energy = -1.0;
  mass   = -1.0;

  // get instance
  static ParticleReco* pRec = ParticleReco::instance();
  float mass    = pRec->invariantMass();
  float energy  = pRec->totalEnergy();

  // compute the momentum "p"
  float p = sqrt( pow(energy,2) - pow(mass,2) );

  // compute distance "d" of the decay point from the origin
  float d = sqrt( pow(ev.getX(),2) + pow(ev.getY(),2) + pow(ev.getZ(),2) );

  // compute the decay proper time
  time = d*mass / (p*Constants::lightVelocity);

  return;

}


float ProperTime::totalEnergy() {
  // check for new event and return result
  check();
  return energy;
}


float ProperTime::invariantMass() {
  // check for new event and return result
  check();
  return mass;
}


ProperTime::ParticleType ProperTime::particleType() {
  // check for new event and return result
  check();
  return type;
}

double ProperTime::decayTime() {
  // check for new event and return result
  check();
  return time;
}
