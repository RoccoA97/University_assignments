#ifndef MassMean_h
#define MassMean_h

class Event;

class MassMean {

	public:

  	MassMean( float min, float max );					// mass range
  	~MassMean();

		// add data from a new event
  	bool add( const Event& ev );
		// compute mean and rms
  	void compute();

		// return number of accepted events
  	unsigned int getAcceptedEvent() const;
		// return mean mass
  	double getMassMean() const;
		// return rms  mass
  	double getMassRMS() const;

 	private:

  	float minMass; // min mass
  	float maxMass; // max mass

  	unsigned int acceptedEvent;			// number of accepted events
  	double massSum;									// sum of masses
  	double massSquare;							// sum of masses square

  	double massMean;			// mean mass
  	double massRMS;				// rms  mass

};

#endif
