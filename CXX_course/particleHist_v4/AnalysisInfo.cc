#include "AnalysisInfo.h"

using namespace std;

// value to return for keys not found
string AnalysisInfo::defaultString = "";


AnalysisInfo::AnalysisInfo( int argc, char* argv[] ) {
  args.resize( argc );
  int iarg;
  for ( iarg = 1; iarg < argc; ++iarg ) args[iarg] = argv[iarg];
}


AnalysisInfo::~AnalysisInfo() {
}


// get the list of all words
const vector<string>& AnalysisInfo::argList() const {
  return args;
}


// get the word coming after the word 'key'
const string& AnalysisInfo::value( const string& key ) const {
  // loop over words
  int i = 0;
  int n = args.size();
  while ( i < n ) {
    // if word is equal to key, return next word
    if ( args[i++] == key ) return args[i];
  }
  // if key not found, return a default string
  return defaultString;
}


// check if the word 'key' is present
bool AnalysisInfo::contains( const string& key ) const {
  // loop over words
  int i = 0;
  int n = args.size();
  while ( i < n ) {
    // if word is equal to key, return true
    if ( args[i++] == key ) return true;
  }
  // if key not found, return false
  return false;
}

