#ifndef ElementReco_h
#define ElementReco_h

#include "../AnalysisFramework/AnalysisSteering.h"
#include "../AnalysisObjects/ProperTime.h"
#include "../AnalysisObjects/LifetimeFit.h"
#include "../util/include/ActiveObserver.h"

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
    virtual void   endJob();
    // function to be called for each event
    virtual void update( const Event& ev );

  private:

    // vector to handle time values
    std::vector<Particle*> timeVector;

    ParticleLifetime           ( const ParticleLifetime& x );
    ParticleLifetime& operator=( const ParticleLifetime& x );

    void pCreate( const std::string& name,
                  double minMass, double maxMass,
                  double minTime, double maxTime,
                  double minScan, double maxScan, double scanStep );

};

#endif
