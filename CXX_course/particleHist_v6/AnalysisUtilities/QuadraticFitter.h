#ifndef QuadraticFitter_h
#define QuadraticFitter_h

class QuadraticFitter {

  // class to fit y = a+bx+cx^2

 public:

  QuadraticFitter();
  virtual ~QuadraticFitter();

  void add( double x, double y );

  double a();
  double b();
  double c();

  void clear();
  void reset() const;

  void update() const;

 private:

  // dummy copy constructor and assignment to prevent unadvertent copy
  QuadraticFitter           ( const QuadraticFitter& x );
  QuadraticFitter& operator=( const QuadraticFitter& x );

  double sumx0y0;
  double sumx1y0;
  double sumx2y0;
  double sumx3y0;
  double sumx4y0;
  double sumx0y1;
  double sumx1y1;
  double sumx2y1;

  mutable bool outdated;
  mutable double ac;
  mutable double bc;
  mutable double cc;

};

#endif

