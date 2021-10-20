#ifndef ParticleReco_h
#define ParticleReco_h

#include "Event.h"

#include "util/include/Singleton.h"
#include "util/include/LazyObserver.h"

class ParticleReco: public Singleton<ParticleReco>,
                    public LazyObserver<Event> {

  friend class Singleton<ParticleReco>;

 public:

  // particle types
  enum ParticleType { K0, Lambda0, unknown };

  // recompute informations for new event
  virtual void update( const Event& ev );

  // return particle energy
  float totalEnergy();
  // return particle mass
  float invariantMass();
  // return particle type
  ParticleType particleType();

 private:

  // private constructor and destructor for singleton
  ParticleReco();
  ~ParticleReco();

  // dummy copy constructor and assignment to prevent unadvertent copy
  ParticleReco           ( const ParticleReco& x );
  ParticleReco& operator=( const ParticleReco& x );

  // particle type
  ParticleType type;
  // particle energy
  float energy;
  // particle mass
  float mass;

};

#endif

