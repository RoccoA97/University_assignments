#include "QuadraticFitter.h"

QuadraticFitter::QuadraticFitter():
 sumx0y0( 0.0 ),
 sumx1y0( 0.0 ),
 sumx2y0( 0.0 ),
 sumx3y0( 0.0 ),
 sumx4y0( 0.0 ),
 sumx0y1( 0.0 ),
 sumx1y1( 0.0 ),
 sumx2y1( 0.0 ) {
  reset();
}


QuadraticFitter::~QuadraticFitter() {
}

void QuadraticFitter::add( double x, double y ) {
  outdated = true;
  double x2 = x  * x;
  double x3 = x2 * x;
  double x4 = x3 * x;
  sumx0y0 += 1.0;
  sumx1y0 += x;
  sumx2y0 += x2;
  sumx3y0 += x3;
  sumx4y0 += x4;
  sumx0y1 += y;
  sumx1y1 += y * x;
  sumx2y1 += y * x2;
  return;
}

double QuadraticFitter::a() {
  if ( outdated ) update();
  return ac;
}

double QuadraticFitter::b() {
  if ( outdated ) update();
  return bc;
}

double QuadraticFitter::c() {
  if ( outdated ) update();
  return cc;
}

void QuadraticFitter::clear() {
  sumx0y0 = 0.0;
  sumx1y0 = 0.0;
  sumx2y0 = 0.0;
  sumx3y0 = 0.0;
  sumx4y0 = 0.0;
  sumx0y1 = 0.0;
  sumx1y1 = 0.0;
  sumx2y1 = 0.0;
  reset();
  return;
}

void QuadraticFitter::reset() const {
  outdated = true;
  ac = bc = cc = 0.0;
  return;
}

void QuadraticFitter::update() const {
  double disc = ( sumx0y0 * sumx2y0 * sumx4y0 ) +
                ( sumx1y0 * sumx3y0 * sumx2y0 ) +
                ( sumx2y0 * sumx1y0 * sumx3y0 ) -
                ( sumx0y0 * sumx3y0 * sumx3y0 ) -
                ( sumx1y0 * sumx1y0 * sumx4y0 ) -
                ( sumx2y0 * sumx2y0 * sumx2y0 );
   double ca0 = ( ( sumx2y0 * sumx4y0 ) - ( sumx3y0 * sumx3y0 ) ) / disc;
   double ca1 = ( ( sumx3y0 * sumx2y0 ) - ( sumx1y0 * sumx4y0 ) ) / disc;
   double ca2 = ( ( sumx1y0 * sumx3y0 ) - ( sumx2y0 * sumx2y0 ) ) / disc;
   double cb0 = ( ( sumx2y0 * sumx3y0 ) - ( sumx1y0 * sumx4y0 ) ) / disc;
   double cb1 = ( ( sumx0y0 * sumx4y0 ) - ( sumx2y0 * sumx2y0 ) ) / disc;
   double cb2 = ( ( sumx1y0 * sumx2y0 ) - ( sumx0y0 * sumx3y0 ) ) / disc;
   double cc0 = ( ( sumx1y0 * sumx3y0 ) - ( sumx2y0 * sumx2y0 ) ) / disc;
   double cc1 = ( ( sumx2y0 * sumx1y0 ) - ( sumx0y0 * sumx3y0 ) ) / disc;
   double cc2 = ( ( sumx0y0 * sumx2y0 ) - ( sumx1y0 * sumx1y0 ) ) / disc;
   ac = ( ca0 * sumx0y1 ) + ( ca1 * sumx1y1 ) + ( ca2 * sumx2y1 );
   bc = ( cb0 * sumx0y1 ) + ( cb1 * sumx1y1 ) + ( cb2 * sumx2y1 );
   cc = ( cc0 * sumx0y1 ) + ( cc1 * sumx1y1 ) + ( cc2 * sumx2y1 );
   outdated = false;
}

