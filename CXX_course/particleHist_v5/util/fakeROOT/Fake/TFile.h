#ifndef TFile_h
#define TFile_h

#include "TDirectory.h"
#include <iostream>
#include <string>

class TH1F;

class TFile: public TDirectory {

 public:

  TFile( const char* name, const char* mode = "READ" );
  virtual ~TFile();

  void Close();

 private:

  static std::ios* open( const char* name, const std::string& mode );
  static bool check    ( const char* name, const std::string& mode );
  virtual void update();

  TFile           ( const TFile& x );
  TFile& operator=( const TFile& x );

};

#endif

