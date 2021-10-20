#include "TH1F.h"
#include "TDirectory.h"
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

TH1F::TH1F( const char* name, const char* title,
            int nbin, float xmin, float xmax ):
 FakeTH1F( name, title, nbin, xmin, xmax ) {
  gDirectory->hmap[hname] = this;
}


TH1F::TH1F( FakeTH1F* fhp ):
 FakeTH1F( *fhp ) {
  gDirectory->hmap[hname] = this;
}


TH1F::~TH1F() {
}


void TH1F::Fill( float x, float w ) {
  int i = ( x + bwid - hmin ) / bwid;
  if ( i < 1    ) i = 0;
  if ( i > hbin ) i = hbin + 1;
  cont[i] += w;
  erro[i] += w * w;
  return;
}


void TH1F::SetBinContent( int i, float x ) {
  if ( i < 1    ) i = 0;
  if ( i > hbin ) i = hbin + 1;
  cont[i] = x;
  return;
}


void TH1F::SetBinError( int i, float e ) {
  if ( i < 1    ) i = 0;
  if ( i > hbin ) i = hbin + 1;
  erro[i] = e * e;
  return;
}


void TH1F::Write() const {
  if ( hname == "" ) return;
  ios* fptr = gDirectory->file();
  ofstream* file = dynamic_cast<ofstream*>( fptr );
  if ( ( fptr != 0 ) && ( file == 0 ) ) {
    cout << "file not open for write" << endl;
    return;
  }
  write( file );
  return;
}


const char* TH1F::GetName() const {
  return hname.c_str();
}


int TH1F::GetNbinsX() const {
  return hbin;
}


int TH1F::FindBin( float x ) const {
  if ( x < hmin ) return 0;
  if ( x > hmax ) return hbin + 1;
  return ( hbin * ( x - hmin ) / ( hmax - hmin ) ) + 1;
}


float TH1F::GetBinContent( int i ) const {
  if ( i < 0 ) return 0.0;
  if ( i > hbin + 1 ) return 0.0;
  return cont[i];
}


float TH1F::GetBinError( int i ) const {
  if ( i < 0 ) return 0.0;
  if ( i > hbin + 1 ) return 0.0;
  return sqrt( erro[i] );
}


TH1F* TH1F::Clone() const {
  if ( hname == "" ) return 0;
  TH1F* h = new TH1F( hname.c_str(), htitle.c_str(), hbin, hmin, hmax );
  int n = hbin + 2;
  int i;
  for ( i = 0; i < n; ++i ) {
    h->cont[i] = cont[i];
    h->erro[i] = erro[i];
  }
  return h;
}


TH1F* TH1F::get( std::ifstream* is ) {
  FakeTH1F* fhp = FakeTH1F::get( is );
  if ( fhp == 0 ) return 0;
  return new TH1F( fhp );
}

