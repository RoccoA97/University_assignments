#include <iostream>
#include <string>
using namespace std;

int main( int argc, char* argv[] ) {

  string* args = new string[argc];
  int iarg;
  for ( iarg = 0; iarg < argc; ++iarg ) args[iarg] = argv[iarg];

  string path = args[1];

  bool found = false;
  for ( iarg = 0; iarg < argc; ++iarg ) {
    if ( args[iarg] == "--cflags" ) {
      if ( found ) cout << " ";
      cout << "-I " << path << "/Core -I " << path << "/Fake";
      found = true;
    }
    if ( args[iarg] == "--libs" ) {
      if ( found ) cout << " ";
      cout << "-L " << path << "/lib -lFakeROOT -lFakeCore";
      found = true;
    }
  }
  cout << endl;

  return 0;

}

