#ifndef ElementReco_h
#define ElementReco_h

#include "AnalysisSteering.h"
#include "ProperTime.h"
#include "LifetimeFit.h"

#include "util/include/ActiveObserver.h"

#include <string>
#include <vector>

class TH1F;
class Event;
class MassMean;

class ParticleLifetime: public AnalysisSteering,
                        public ActiveObserver<Event> {

  public:

    ParticleLifetime( const AnalysisInfo* info );
    virtual ~ParticleLifetime();

    struct Particle {
      std::string name;
      LifetimeFit* lifetimeFitPtr;
      TH1F* rootHist;
    };

    // function to be called at execution start
    virtual void beginJob();
    // function to be called at execution end
    virtual void endJob();
    // function to be called for each event
    virtual void update( const Event& ev );

  private:

    // vector to handle time values
    std::vector<Particle*> timeVector;

    ParticleLifetime           ( const ParticleLifetime& x );
    ParticleLifetime& operator=( const ParticleLifetime& x );

    void pCreate( const std::string& name, float min, float max, float timeMin, float timeMax );

};

#endif
