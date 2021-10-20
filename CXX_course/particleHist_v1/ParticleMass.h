#ifndef ElementReco_h
#define ElementReco_h

#include <string>
#include <vector>

#include "AnalysisSteering.h"

class TH1F;
class Event;
class MassMean;


class ParticleMass: public AnalysisSteering {

  public:

    ParticleMass();
    virtual ~ParticleMass();

    struct Particle {
      std::string name;
      MassMean* massMeanPtr;
      TH1F* rootHist;
    };

    // function to be called at execution start
    virtual void beginJob();
    // function to be called at execution end
    virtual void   endJob();
    // function to be called for each event
    virtual void process( const Event& ev );

  private:

    // set of Particle pointers
    std::vector<Particle*> massVector;

    ParticleMass           ( const ParticleMass& x );
    ParticleMass& operator=( const ParticleMass& x );

    void pCreate( const std::string& name, float min, float max );

};

#endif
