#ifndef AnalysisSteering_h
#define AnalysisSteering_h

class Event;

class AnalysisSteering {

  public:

    AnalysisSteering();
    virtual ~AnalysisSteering();

    // function to be called at execution start
    virtual void beginJob() = 0;
    // function to be called at execution end
    virtual void   endJob() = 0;
    // function to be called for each event
    virtual void process( const Event& ev ) = 0;

  private:

    AnalysisSteering           ( const AnalysisSteering& x );
    AnalysisSteering& operator=( const AnalysisSteering& x );

};

#endif
