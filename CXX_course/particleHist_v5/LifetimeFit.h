#ifndef LifetimeFit_h
#define LifetimeFit_h

class Event;

class LifetimeFit {

	public:

  		LifetimeFit( 	float min, float max,
										double minTime, double maxTime,
										double minScan, double maxScan, double scanStep );
  		~LifetimeFit();

			// add data from a new event
  		bool add( const Event& ev );
			// compute mean and rms
  		void compute();

			// return number of accepted events
  		unsigned int nEvents() const;

			// return the mean life time and error
			double lifeTime() const;
			double lifeTimeError() const;



 	private:

  		float minMass; // min mass
  		float maxMass; // max mass

			double minTime;	// min time
			double maxTime;	// max time
			double minScan;
			double maxScan;
			double scanStep;

			std::vector<double> decayTimes;

			double ptLifeTimeMean;
			double ptLifeTimeError;


      // dummy copy constructor and assignment to prevent unadvertent copy
      LifetimeFit           ( const LifetimeFit& x );
      LifetimeFit& operator=( const LifetimeFit& x );

};

#endif
