#ifndef EventSim_h
#define EventSim_h

#include "EventSource.h"

#include <math.h>
#include <stdlib.h>

class Event;

class EventSim: public EventSource {

 public:

  // simulate n events with random seed s
  EventSim( unsigned int n, unsigned int s );
  virtual ~EventSim();

 private:

  // number of events to simulate
  int nMax;

  // last event identifier
  int evId;

  // get an event
  virtual const Event* get();

  // generate an event
  const Event* generate();

  // constants
  static const double massPion;          // masses (GeV/c^2)
  static const double massK0;
  static const double massLambda0;
  static const double massProton;
  static const double lightVelocity;     // c (cm/ps)
  static const double lifeK0;            // lifetimes (ps)
  static const double lifeLambda0;
  static const float  fractionK0;        // produced fractions
  static const float  fractionLambda0;
  static const int    maxParticlesBG;    // max particles in background
//  static const int    maxParticlesMix;   // max particles in background mix
  static const float  beamDivK0;         // beam divergence
  static const float  beamDivLambda0;
  static const float  partDivBG;         // background divergence
  static const float  maxFactK0;         // max cos_complement factor
  static const float  maxFactLambda0;
  static const float  maxFactBG;
  static const float  meanEnergyK0;      // mean particle energy
  static const float  meanEnergyLambda0;
  static const float  meanEnergyBG;
  static const float   rmsEnergyK0;      // particle energy dispersion
  static const float   rmsEnergyLambda0;
  static const float   rmsEnergyBG;
  static const float   minEnergyK0;      // minimum particle energy
  static const float   minEnergyLambda0;
  static const float   minEnergyBG;

    // K0 Lambda0
  
  // generate a K0
  const Event* genK0     ( int id );
  // generate a Lambda0
  const Event* genLambda0( int id );
  // generate background
  const Event* genBG     ( int id );

  // generate decay
  const Event* genDecay( int id,
                         double mass, double life, double energy,
                         double bCos, double bPhi,
                         double mPos, double mNeg );

  // add background particle
  void addBG( Event* ev );

  // random number generation

  // uniform 
  static double randF() { return         random() *   1.0         /
                                                      RAND_MAX     ; }
  static double randF( float max ) {
                          return         random() *   max         /
                                                      RAND_MAX     ; }
  static double randF( float min, float max ) {
                          return min + ( random() * ( max - min ) /
                                                      RAND_MAX    ); }

  // exponential
  static double randE()  { return     -log( randF() ); }
  static double randE( float c )  {
                           return -c * log( randF() ); }

  // gaussian
  static double randG() { return                  sqrt( -log( randF() ) ) *
                                                         sin( randPhi() )    ; }
  static double randG( float sigma ) {
                          return          sigma * sqrt( -log( randF() ) ) *
                                                         sin( randPhi() )    ; }
  static double randG( float mean, float sigma ) {
                          return mean + ( sigma * sqrt( -log( randF() ) ) *
                                                         sin( randPhi() ) )  ; }

  // positive gaussian
  static double randH() { return                  sqrt( -log( randF() ) ) *
                                                         sin( randF( M_PI ) ); }

  // energy
  static double randE( float mean, float sigma, float min ) {
    double e = -1.0;
    while ( e < min ) e = randG( mean, sigma );
    return e;
  }

  // phi angle
  static double randPhi() { return randF( 2.0 * M_PI ); }

  // cosine
  static double randCos( float cdiv, float fact ) {
    double r = 0.0;
    double c = cos( cdiv );
    double m = complement( c, fact );
    while ( r < m ) r = complement( c, randH() );
    return r;
  }

  // charge
  static int randC() { return ( randF() > 0.5 ? +1 : -1 ); }

  // cosine complement
  static float complement( float c, float f ) {
    return 1.0 - ( ( 1.0 - c ) * f );
  }

  // cartesian 3-vector for positions / momenta
  class Vector3 {

   public:

    // create point in the origin / null vector
    Vector3():
     xv( 0.0 ),
     yv( 0.0 ),
     zv( 0.0 ),
     ms( -1.0 ),
     mv( -1.0 ) {}

    // create point / momentum with coordinates / components "x,y,z"
    Vector3( double x, double y, double z ):
     xv( x ),
     yv( y ),
     zv( z ),
     ms( -1.0 ),
     mv( -1.0 ) {}

    // copy constructor and assignment
    Vector3           ( const Vector3& v ) { copy( v ); }
    Vector3& operator=( const Vector3& v ) { copy( v ); return *this; }

    // destructor
    virtual ~Vector3() {}

    // return coordinates / components
    double x() const { return xv; }
    double y() const { return yv; }
    double z() const { return zv; }

    // return ( squared ) modulus
    double msq() const { if ( ms < 0.0 ) update(); return ms; }
    double mod() const { if ( ms < 0.0 ) update(); return mv; }

    // scale vector to modulus "m"
    Vector3 norm( double m = 1.0 ) const {
      double f = m / mod();
      return Vector3( xv * f, yv * f, zv * f );
    }

    // sum vectors
    Vector3 operator+( const Vector3& v ) const {
      return Vector3( xv + v.xv, yv + v.yv, zv + v.zv );
    }

    // subtract vector
    Vector3 operator-( const Vector3& v ) const {
      return Vector3( xv - v.xv, yv - v.yv, zv - v.zv );
    }

    // reverse vector
    Vector3 operator-() const {
      return Vector3( -xv, -yv, -zv );
    }

    // scale vector by factor "f"
    Vector3 operator*( double f ) const {
      return Vector3( f * xv, f * yv, f * zv );
    }

    // scale vector by factor "1/f"
    Vector3 operator/( double f ) const {
      return *this * ( 1.0 / f );
    }

    // scalar product
    double operator|( const Vector3& v ) const {
      return ( xv * v.xv ) + ( yv * v.yv ) + ( zv * v.zv );
    }

    // component along vector "v"
    double operator%( const Vector3& v ) const {
      return ( *this | v ) / v.mod();
    }

    // vector product
    Vector3 operator^( const Vector3& v ) const {
      return Vector3( ( yv * v.zv ) - ( zv * v.yv ),
                      ( zv * v.xv ) - ( xv * v.zv ),
                      ( xv * v.yv ) - ( yv * v.xv ) );
    }

    // vector component along vector "v"
    Vector3 operator||( const Vector3& v ) const {
      return v * ( ( v | *this ) / v.msq() );
    }

    // vector component normal to vector "v"
    Vector3 operator&&( const Vector3& v ) const {
      return *this - ( *this || v );
    }

    // squared modulus of vector with components "x,y.z"
    static double msq( double x, double y, double z ) {
      return ( x * x ) + ( y * y ) + ( z * z );
    }

    // modulus of vector with components "x,y.z"
    static double mod( double x, double y, double z ) {
      return sqrt( msq( x, y, z ) );
    }

   private:

    // coordinates / components
    double xv;
    double yv;
    double zv;
    mutable double ms;
    mutable double mv;

    // copy vector
    void copy( const Vector3& v ) {
      xv = v.xv;
      yv = v.yv;
      zv = v.zv;
      ms = -1.0;
      return;
    }

    // update ( squared ) modulus
    void update() const {
      ms = msq( xv, yv, zv );
      mv = sqrt( ms );
      return;
    }

  };

  // Lorentz 4-vector
  class Vector4 {

   public:

    // create null vector
    Vector4():
    pm( 0.0 ) { ee = 0.0; }

    // create vector with 3-momentum "px,py,pz" and mass "im"
    Vector4( double px, double py, double pz, double im ):
     vp( px, py, pz ),
     pm( im ) { ee = -1.0; }

    // create vector with 3-momentum "pv" and mass "im"
    Vector4( const Vector3& pv, double im ):
     vp( pv ),
     pm( im ) { ee = -1.0; }

    // copy constructor and assignment
    Vector4           ( const Vector4& v ) { copy( v ); }
    Vector4& operator=( const Vector4& v ) { copy( v ); return *this; }

    // destructor
    virtual ~Vector4() {}

    // 3-momentum, energy, mass, velocity (beta) and boost (gamma)
    const Vector3& p() const { return vp; }
    double energy() const { if ( ee < 0.0 ) update(); return ee; }
    double beta() const {
      return vp.mod() / energy();
    }
    double gamma() const {
      return energy() / pm;
    }

    // boosted 4-vector along 4-vector "v"
    Vector4 boost( const Vector4& v ) {
      Vector3 u = vp || v.p();
      Vector3 t = vp - u;
      return Vector4( v.p().norm( ( ( vp % v.p() ) + ( v.beta() * energy() ) )
                      * v.gamma() )
                      + t, pm );
    }

    // ( squared ) energy of 4-vector with components "px,py,pz" and mass "im"
    static double enersq( double px, double py, double pz, double im ) {
      return ( px * px ) + ( py * py ) + ( pz * pz ) + ( im * im );
    }
    static double energy( double px, double py, double pz, double im ) {
      return sqrt( enersq( px, py, pz, im ) );
    }

    // ( squared ) energy of 4-vector with 3-momentum "pv" and mass "im"
    static double enersq( const Vector3& pv, double im ) {
      return pv.msq() + ( im * im );
    }
    static double energy( const Vector3& pv, double im ) {
      return sqrt( enersq( pv, im ) );
    }

   private:

    // 3-momentum
    Vector3 vp;
    // mass
    double pm;
    // energy
    mutable double ee;

    void copy( const Vector4& v ) {
      vp = v.vp;
      pm = v.pm;
      ee = v.ee;
    }

    void update() const {
      ee = energy( vp, pm );
      return;
    }

  };

  // dummy copy constructor and assignment to prevent unadvertent copy
  EventSim           ( const EventSim& x );
  EventSim& operator=( const EventSim& x );

};

#endif

