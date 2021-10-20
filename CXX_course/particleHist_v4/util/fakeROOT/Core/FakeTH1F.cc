#include "FakeTH1F.h"
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

FakeTH1F::FakeTH1F( const char* name, const char* title,
                    int nbin, float xmin, float xmax ):
 hname ( name ),
 htitle( title ),
 hbin( nbin ),
 hmin( xmin ),
 hmax( xmax ) {
  int nval = nbin + 2;
  cont = new float[nval];
  erro = new float[nval];
  while ( nval-- ) cont[nval] = erro[nval] = 0.0;
  bwid = ( hmax - hmin ) / hbin;
}


FakeTH1F::FakeTH1F() {
}


FakeTH1F::FakeTH1F( FakeTH1F& fh ):
 hname ( fh.hname ),
 htitle( fh.htitle ),
 hbin( fh.hbin ),
 hmin( fh.hmin ),
 hmax( fh.hmax ),
 cont ( fh.cont ),
 erro ( fh.erro ) {
  fh.cont = 0;
  fh.erro = 0;
}


FakeTH1F::~FakeTH1F() {
  delete[] cont;
  delete[] erro;
}


void FakeTH1F::write( std::ofstream* os ) const {
  write( hname , os );
  write( htitle, os );
  write( hbin, hmin, hmax, cont, erro, os );
  return;
}


FakeTH1F* FakeTH1F::get( std::ifstream* is ) {
  string name  = read( is );
  if ( name == "" ) return 0;
  string title = read( is );
  int n;
  float a;
  float b;
  is->read( reinterpret_cast<char*>(      &n ),     sizeof(        n ) );
  is->read( reinterpret_cast<char*>(      &a ),     sizeof(        a ) );
  is->read( reinterpret_cast<char*>(      &b ),     sizeof(        b ) );
  FakeTH1F* h = new FakeTH1F( name.c_str(), title.c_str(), n, a, b );
  n += 2;
  is->read( reinterpret_cast<char*>( h->cont ), n * sizeof( *h->cont ) );
  is->read( reinterpret_cast<char*>( h->erro ), n * sizeof( *h->erro ) );
  return h;
}


void FakeTH1F::write( const std::string& s, std::ofstream* os ) {
  if ( os == 0 ) {
    cout << s << endl;
    return;
  }
  int n = s.length();
  os->write( reinterpret_cast<char*>( &n ), sizeof( n ) );
  os->write( s.c_str(), n );
  return;
}


void FakeTH1F::write( int n, float a, float b,
                      float* c, float* e, std::ofstream* os ) {
  int i;
  if ( os == 0 ) {
    cout << n << " " << a << " " << b << endl;
    for ( i = 0; i < n; ++i ) cout << c[i] << " " << sqrt( e[i] ) << endl;
    return;
  }
  os->write( reinterpret_cast<char*>( &n ),     sizeof(  n ) );
  os->write( reinterpret_cast<char*>( &a ),     sizeof(  a ) );
  os->write( reinterpret_cast<char*>( &b ),     sizeof(  b ) );
  n += 2;
  os->write( reinterpret_cast<char*>(  c ), n * sizeof( *c ) );
  os->write( reinterpret_cast<char*>(  e ), n * sizeof( *e ) );
  return;
}


string FakeTH1F::read( std::ifstream* is ) {
  int l;
  char* p;
  if ( is->read( reinterpret_cast<char*>( &l ), sizeof( l ) ) )
       p = new char[l + 1];
  else return "";
  is->read( p, l );
  p[l] = '\0';
  string s( p );
  delete p;
  return s;
}

