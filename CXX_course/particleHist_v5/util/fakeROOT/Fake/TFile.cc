#include "TFile.h"
#include "TH1F.h"

using namespace std;

TFile::TFile( const char* name, const char* mode ) {
  string ms = mode;
  fp = open( name, ms );
  cd();
}


TFile::~TFile() {
  Close();
}


void TFile::Close() {
  if ( fp == &cin ) fp = 0;
  if ( fp != 0 ) delete fp;
  fp = 0;
}


ios* TFile::open( const char* name, const string& mode ) {
  if ( check ( name, mode ) ) return &cin;
  ios* fp;
  if ( ( mode ==      "NEW" ) ||
       ( mode ==   "CREATE" ) ||
       ( mode == "RECREATE" ) ) fp = new std::ofstream( name, ios::binary );
  else
  if (   mode ==   "UPDATE"   ) fp = new std::ofstream( name, ios::binary |
                                         std::ofstream::in  |
                                         std::ofstream::out |
                                         std::ofstream::app );
  else                          fp = new std::ifstream( name, ios::binary );
  return fp;
}


bool TFile::check( const char* name, const std::string& mode ) {
  if ( ( mode == "RECREATE" ) ||
       ( mode ==   "UPDATE" ) ) return false;
  ifstream file( name );
  bool status = static_cast<bool>( file );
  if ( ( mode ==      "NEW" ) ||
       ( mode ==   "CREATE" ) ) return status;
  return !status;
}


void TFile::update() {
  ifstream* file = dynamic_cast<ifstream*>( fp );
  if ( file == 0 ) return;
  TDirectory* currentDir = gDirectory;
  cd();
  while ( TH1F::get( file ) );
  currentDir->cd();
  return;
}

