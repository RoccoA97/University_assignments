#ifndef Particle_h
#define Particle_h

struct Particle {
	int charge;
	float px;
	float py;
	float pz;
	};
	
struct Event {
	int eventId;
	int nPart;
	float x;
	float y;
	float z;
	Particle** pt;
	};
	
#endif
