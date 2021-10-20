#ifndef FakeTH1F_h
#define FakeTH1F_h

#include <iostream>
#include <fstream>
#include <string>

class FakeTH1F {

  friend class TH1F;

 public:

  FakeTH1F( const char* name, const char* title,
        int nbin, float xmin, float xmax );
  virtual ~FakeTH1F();

  static FakeTH1F* get( std::ifstream* is );

 protected:

  std::string hname;
  std::string htitle;

  int   hbin;
  float hmin;
  float hmax;
  float bwid;
  float* cont;
  float* erro;

  void write( std::ofstream* os ) const;

  FakeTH1F();
  FakeTH1F( FakeTH1F& fh );

 private:

  static void write( const std::string& s, std::ofstream* os );
  static void write( int n, float a, float b,
                     float* c, float* e, std::ofstream* os );
  static std::string read( std::ifstream* is );

  FakeTH1F           ( const FakeTH1F& x );
  FakeTH1F& operator=( const FakeTH1F& x );

};

#endif

