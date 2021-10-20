#include "Random.h"

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

int main( int argc, char* argv[] ) {

  Random::verbosity = 3;

  string strn = argv[1];
  string prob = argv[2];
  string str1 = argv[3];
  string str2 = argv[4];

  stringstream sstr;
  sstr.clear();
  sstr.str( str1 );
  float par1;
  sstr >> par1;
  sstr.clear();
  sstr.str( str2 );
  float par2;
  sstr >> par2;
  sstr.clear();
  sstr.str( strn );
  int n;
  sstr >> n;

  switch ( *prob.c_str() ) {
  case 'f':
    while ( n-- ) cout << Random::flat ( par1, par2 ) << endl;
    break;
  case 'g':
    while ( n-- ) cout << Random::gauss( par1, par2 ) << endl;
    break;
  }

  return 0;

}

