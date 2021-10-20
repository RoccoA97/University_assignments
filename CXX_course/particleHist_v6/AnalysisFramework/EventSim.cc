#include "EventSim.h"
#include "Event.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

// constants
const double EventSim::massPion          =   0.1395706; // masses (GeV/c^2)
const double EventSim::massK0            =   0.497611;
const double EventSim::massLambda0       =   1.115683;
const double EventSim::massProton        =   0.938272;
const double EventSim::lightVelocity     =   0.029979246; // c (cm/ps)
const double EventSim::lifeK0            =  89.54; // lifetimes (ps)
const double EventSim::lifeLambda0       = 263.2 ;
const float  EventSim::fractionK0        =   0.5 ; // produced fraction
const float  EventSim::fractionLambda0   =   0.3 + EventSim::fractionK0;
const int    EventSim::maxParticlesBG    =   5   ; // max particles in bg
//const int    EventSim::maxParticlesMix   =   1   ; // max particles in bg mix
const float  EventSim::beamDivK0         =   0.1 ; // beam divergence
const float  EventSim::beamDivLambda0    =   0.1 ;
const float  EventSim::partDivBG         =   0.1 ; // background divergence
const float  EventSim::maxFactK0         =   2.5 ; // max cos_complement factor
const float  EventSim::maxFactLambda0    =   2.5 ;
const float  EventSim::maxFactBG         =   2.5 ; // max cos_complement factor
const float  EventSim::meanEnergyK0      =   5.0 ; // mean particle energy
const float  EventSim::meanEnergyLambda0 =   5.0 ;
const float  EventSim::meanEnergyBG      =   1.0 ;
const float  EventSim:: rmsEnergyK0      =   2.0 ; // particle energy dispersion
const float  EventSim:: rmsEnergyLambda0 =   2.0 ;
const float  EventSim:: rmsEnergyBG      =   0.5 ;
const float  EventSim:: minEnergyK0      =   2.0 ; // minimum particle energy
const float  EventSim:: minEnergyLambda0 =   2.0 ;
const float  EventSim:: minEnergyBG      =   0.2 ;

// simulate data with random seed
EventSim::EventSim( unsigned int n, unsigned int s ):
 nMax( n ),
 evId( 0 ) {
  srandom( s );
}


EventSim::~EventSim() {
}


// get an event
const Event* EventSim::get() {
  if ( nMax-- ) return generate();
  return 0;
}


// generate an event
const Event* EventSim::generate() {

  // set event id
  evId += ceil( randE( 10.0 ) );

  // choose channel
  float frac = randF();
  if ( frac < fractionK0 )
  return genK0     ( evId );
  if ( frac < fractionLambda0 )
  return genLambda0( evId );
  return genBG     ( evId );

}


// generate a K0
const Event* EventSim::genK0( int id ) {

  // direction
  double cosBeam = randCos( beamDivK0, maxFactK0 );
  double phiBeam = randPhi();

  // energy
  double energy = randE( meanEnergyK0,
                          rmsEnergyK0,
                          minEnergyK0 );

  // decay
  return genDecay( id, massK0, lifeK0,
                   energy, cosBeam, phiBeam,
                   massPion, massPion );

}


// generate a Lambda0
const Event* EventSim::genLambda0( int id ) {

  // direction
  double cosBeam = randCos( beamDivLambda0, maxFactLambda0 );
  double phiBeam = randPhi();

  // energy
  double energy = randE( meanEnergyLambda0,
                          rmsEnergyLambda0,
                          minEnergyLambda0 );

  // decay
  return genDecay( id, massLambda0, lifeLambda0,
                   energy, cosBeam, phiBeam,
                   massProton, massPion );

}


// generate background
const Event* EventSim::genBG( int id ) {

  // vertex point
  static float svx = 0.2;
  static float svy = 0.2;
  static float svz = 0.2;
  Event* ev = new Event( id, randG( svx ), randG( svy ), randG( svz ) );

  // number of particles
  int n = ceil( randE() );
  if ( n > maxParticlesBG ) n = maxParticlesBG;
  while ( n-- ) addBG( ev );

  return ev;

}


// generate decay
const Event* EventSim::genDecay( int id,
                                 double mass, double life, double energy,
                                 double bCos, double bPhi,
                                 double mPos, double mNeg ) {

  // beam direction
  double bSin = sqrt( 1.0 - ( bCos * bCos ) );

  // boost and 4-momentum
  double gamma = energy / mass;
  double pb = sqrt( ( energy * energy ) - ( mass * mass ) );
  Vector4 p0( pb * bSin * cos( bPhi ),
              pb * bSin * sin( bPhi ),
              pb * bCos, mass );

  // flight distance and decay point
  double dist = randE( gamma * lightVelocity * life );
  Vector3 pv = p0.p().norm();

  // create a new event
  static float pRMS = 0.05;
  double xVtx = randG( dist * pv.x(), pRMS );
  double yVtx = randG( dist * pv.y(), pRMS );
  double zVtx = randG( dist * pv.z(), pRMS );
  Event* ev = new Event( id, xVtx, yVtx, zVtx );

  // decay in rest frame
  double dCos = randF( -1.0, 1.0 );
  double dSin = sqrt( 1.0 - ( dCos * dCos ) );
  double dPhi = randPhi();
  double pd = sqrt( ( mass * mass / 4.0 ) +
                    ( ( pow( mPos, 4 ) + pow( mNeg, 4 ) -
                        ( ( 2 * mPos * mPos * mNeg * mNeg ) +
                          ( 2 * mass * mass * mPos * mPos ) +
                          ( 2 * mass * mass * mNeg * mNeg ) ) ) /
                     ( mass * mass * 4.0 ) ) );
  double px = pd * dSin * cos( dPhi );
  double py = pd * dSin * sin( dPhi );
  double pz = pd * dCos;

  // decay in lab. frame
  static float fRMS = 0.01;
  Vector3 pPos = Vector4( randG(  px, px * px * fRMS ),
                          randG(  py, py * py * fRMS ),
                          randG(  pz, pz * pz * fRMS ),
                           mPos ).boost( p0 ).p();
  Vector3 pNeg = Vector4( randG( -px, px * px * fRMS ),
                          randG( -py, py * py * fRMS ),
                          randG( -pz, pz * pz * fRMS ),
                           mNeg ).boost( p0 ).p();
  Vector3* pDaug1;
  Vector3* pDaug2;
  int charge = randC();

  // add particles to event
  if ( charge > 0 ) {
    pDaug1 = &pPos;
    pDaug2 = &pNeg;
  }
  else {
    pDaug1 = &pNeg;
    pDaug2 = &pPos;
  }
  ev->add( pDaug1->x(), pDaug1->y(), pDaug1->z(),  charge );
  ev->add( pDaug2->x(), pDaug2->y(), pDaug2->z(), -charge );

  // add background
//  int n = 0.2 * ceil( randE() );
//  if ( n > maxParticlesMix ) n = maxParticlesMix;
//  while ( n-- ) addBG( ev );

  return ev;

}


// add background particle
void EventSim::addBG( Event* ev ) {
  double cosPart = randCos( partDivBG, maxFactBG );
  double phiPart = randPhi();
  double sinPart = sqrt( 1.0 - ( cosPart * cosPart ) );
  double energy = randE( meanEnergyBG,
                          rmsEnergyBG,
                          minEnergyBG );
  int charge = randC();
  ev->add( energy * sinPart * cos( phiPart ),
           energy * sinPart * sin( phiPart ),
           energy * cosPart,
           charge );
  return;

}
