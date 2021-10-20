#ifndef LifetimeFit_h
#define LifetimeFit_h

class Event;

class LifetimeFit {

	public:

  		LifetimeFit( float min, float max ); // mass range
  		~LifetimeFit();

			// add data from a new event
  		bool add( const Event& ev );
			// compute mean and rms
  		void compute();

			// return number of accepted events
  		unsigned int getAcceptedEvent() const;

 	private:

  		float minMass; // min mass
  		float maxMass; // max mass

			// number of accepted events
			unsigned int acceptedEvent;

      // dummy copy constructor and assignment to prevent unadvertent copy
      LifetimeFit           ( const LifetimeFit& x );
      LifetimeFit& operator=( const LifetimeFit& x );

};

#endif
