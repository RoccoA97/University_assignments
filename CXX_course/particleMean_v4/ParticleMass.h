#ifndef ElementReco_h
#define ElementReco_h

#include <vector>

#include "AnalysisSteering.h"

class Event;
class MassMean;


class ParticleMass: public AnalysisSteering {

  public:

    ParticleMass();
    virtual ~ParticleMass();

    // function to be called at execution start
    virtual void beginJob();
    // function to be called at execution end
    virtual void   endJob();
    // function to be called for each event
    virtual void process( const Event& ev );

  private:

    // set of MassMean pointers for different hypothesis
    std::vector<MassMean*> massVector;

    ParticleMass           ( const ParticleMass& x );
    ParticleMass& operator=( const ParticleMass& x );

};

#endif
