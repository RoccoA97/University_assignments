#ifndef TH1F_h
#define TH1F_h

#include "FakeTH1F.h"
#include <iostream>
#include <fstream>
#include <string>

class TH1F: public FakeTH1F {

 public:

  TH1F( const char* name, const char* title,
        int nbin, float xmin, float xmax );
  virtual ~TH1F();

  void Fill( float x, float w = 1.0 );
  void SetBinContent( int i, float x );
  void SetBinError  ( int i, float e );
  void Write() const;

  const char* GetName() const;
  int GetNbinsX() const;
  int FindBin( float x ) const;
  float GetBinContent( int i ) const;
  float GetBinError( int i ) const;
  TH1F* Clone() const;

  static TH1F* get( std::ifstream* is );

 private:

  TH1F( FakeTH1F* fhp );

  TH1F           ( const TH1F& x );
  TH1F& operator=( const TH1F& x );

};

#endif

