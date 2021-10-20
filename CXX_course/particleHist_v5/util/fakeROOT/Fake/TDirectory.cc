#include "TDirectory.h"
//#include "TH1F.h"
#include <iostream>
#include <math.h>
using namespace std;

TDirectory::TDirectory(): fp( 0 ) {
}


TDirectory::~TDirectory() {
}


void TDirectory::cd() {
  gDirectory = this;
  return;
}


std::ios* TDirectory::file() const {
  return fp;
}


TH1F* TDirectory::Get( const char* name ) {
  if ( !hmap.size() ) update();
  if ( !hmap.size() ) return 0;
  map<string,TH1F*>::iterator iter = hmap.find( name );
  map<string,TH1F*>::iterator iend = hmap.end();
  if ( iter != iend ) return iter->second;
  return 0;
}


std::vector<std::string> TDirectory::histoList() {
  TDirectory* currentDir = gDirectory;
  cd();
  std::vector<std::string> hl;
  if ( !hmap.size() ) update();
  if ( !hmap.size() ) return hl;
  hl.reserve( hmap.size() );
  map<string,TH1F*>::iterator iter = hmap.begin();
  map<string,TH1F*>::iterator iend = hmap.end();
  while ( iter != iend ) hl.push_back( iter++->first );
  currentDir->cd();
  return hl;
}


void TDirectory::update() {
  return;
}


TDirectory* gDirectory = new TDirectory;

