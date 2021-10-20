#ifndef Random_H
#define Random_H

class Random {

 public:

  virtual ~Random();

  // random number generation functions
  static void setSeed( unsigned int seed );
  static float flat ( double min , double max );
  static float gauss( double mean, double rms );

  // verbosity level
  static int verbosity;

 protected:

  // proability distribution
  enum probability { Flat, Gauss };

  // the object can be created only through a derived object
  Random();

  // generic random number generation functions
  // to be implemented by concrete random generators
  virtual void set( unsigned int seed );
  virtual float generate( probability p, float a, float b );

 private:

  // get a pointer to a random generator,
  // eventually to be assigned a concrete random generator
  static Random*  instance();
  static Random*& instRef();

  Random           ( const Random& x );
  Random& operator=( const Random& x );

};

#endif

