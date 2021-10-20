#ifndef MassMean_h
#define MassMean_h

class Event;

class MassMean {

	public:

  		MassMean( float min, float max ); // mass range
  		~MassMean();

  		bool add( const Event& ev );      // add data from a new event
  		void compute();                   // compute mean and rms

  		unsigned int getAcceptedEvent() const;                   // return number of accepted events
  		double getMassMean() const;                               // return mean mass
  		double getMassRMS() const;	                         // return rms  mass

 	private:

  		float minMass; // min mass
  		float maxMass; // max mass

  		unsigned int acceptedEvent; // number of accepted events
  		double massSum; // sum of masses
  		double massSquare; // sum of masses square

  		double massMean; // mean mass
  		double massRMS; // rms  mass

      // dummy copy constructor and assignment to prevent unadvertent copy
      MassMean           ( const MassMean& x );
      MassMean& operator=( const MassMean& x );

};

#endif
