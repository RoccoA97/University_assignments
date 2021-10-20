#ifndef Event_h
#define Event_h

#include <vector>


class Event {

  public:

    // create an event with number "n" and coordinates x, y, z
    Event( int n, float x, float y, float z );
    ~Event();

    // composite object Particle to hold all information for each particle
    // ( x,y,z momentum components and electric charge )
    struct Particle {
      int charge;
      float px;
      float py;
      float pz;
    };

    // add a particle to the event
    void add( float px, float py, float pz, int charge );

    // get event id.
    int eventNumber() const;
    // get decay point coordinates
    float getX() const;
    float getY() const;
    float getZ() const;
    // get number of particles
    unsigned int nParticles() const;
    // get particle
    const Particle* particle( unsigned int i ) const;

  private:

    // dummy copy constructor and assignment to prevent unadvertent copy
    Event           ( const Event& x );
    Event& operator=( const Event& x );

    // event-specific informations:
    int eventId; // event id
    float x;     // decay point
    float y;
    float z;

    // particles: array of pointers
    std::vector<const Particle*> pt;

};

#endif
