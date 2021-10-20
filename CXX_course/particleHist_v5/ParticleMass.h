#ifndef ElementReco_h
#define ElementReco_h

#include "AnalysisSteering.h"

#include "util/include/ActiveObserver.h"

#include <string>
#include <vector>

class TH1F;
class Event;
class MassMean;

class ParticleMass: public AnalysisSteering,
                    public ActiveObserver<Event> {

  public:

    ParticleMass( const AnalysisInfo* info );
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
    virtual void update( const Event& ev );

  private:

    // set of Bragg curves for different total energies
    std::vector<Particle*> massVector;

    ParticleMass           ( const ParticleMass& x );
    ParticleMass& operator=( const ParticleMass& x );

    void pCreate( const std::string& name, float min, float max );

};

#endif
