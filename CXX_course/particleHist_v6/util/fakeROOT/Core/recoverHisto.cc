#include "TFile.h"
#include "TH1.h"
#include "FakeTH1F.cc"

#include <math.h>

using namespace std;

class H1FRecover: public FakeTH1F {

 public:

  H1FRecover( FakeTH1F* fhp );
  TH1F* recover() const;

private:

};

H1FRecover::H1FRecover( FakeTH1F* fhp ):
 FakeTH1F( *fhp ) {
}

TH1F* H1FRecover::recover() const {
  cout << "recover " << hname << endl;
  TH1F* h = new TH1F( hname.c_str(), htitle.c_str(), hbin, hmin, hmax );
  int n = hbin + 2;
  int i;
  float r = 0.0;
  float s = 0.0;
  for ( i = 0; i < n; ++i ) {
    r += fabs( cont[i] - lround( cont[i] ) );
    s += fabs( cont[i] - erro[i] );
  }
  if ( ( r > 1.0 ) || ( s > 1.0 ) ) {
    for ( i = 0; i < n; ++i ) {
      h->SetBinContent( i, cont[i] );
      h->SetBinError  ( i, sqrt( erro[i] ) );
    }
  }
  else {
    float wbin = ( hmax - hmin ) / hbin;
    float xmin = hmin - ( wbin / 2.0 );
    for ( i = 0; i < n; ++i ) {
      int m = lround( cont[i] );
      while ( m-- ) h->Fill( xmin );
      xmin += wbin;
    }
  }
  return h;
}


//int main( int argc, char* argv[] ) {
void recoverHisto( const char* i_name, const char* o_name ) {

  ifstream* i_file = new ifstream( i_name );

  TDirectory* currentDir = gDirectory;
  TFile* o_file = new TFile ( o_name, "RECREATE" );

  FakeTH1F* h;
  while (( h = FakeTH1F::get( i_file ) )) {
    H1FRecover r( h );
    r.recover()->Write();
  }

  currentDir->cd();
  o_file->Close();
  delete o_file;
  delete i_file;

  return;

}

