#include "TFile.h"
#include "TH1.h"

using namespace std;

bool getName( string& name ) {
  cout << "hist > " << flush;
  cin >> name;
  return static_cast<bool>( cin );
}

void dumpList( const vector<string>& list ) {
  int i;
  int n = list.size();
  for ( i = 0; i < n; ++i ) {
    cout << list[i] << endl;
  }
  return;
}

int main( int argc, char* argv[] ) {

  string name;
  if ( argc >= 2 ) name = argv[1];
  if ( name == "" ) {
    cout << "file > " << flush;
    cin >> name;
  }

  TDirectory* currentDir = gDirectory;
  TFile* file = new TFile ( name.c_str() );
  currentDir->cd();

  vector<string> list = file->histoList();
  dumpList( list );

  string hist;
  while ( getName( hist ) ) {
    if ( hist == "." ) {
      dumpList( list );
      continue;
    }
    TH1F* h = dynamic_cast<TH1F*>( file->Get( hist.c_str() ) );
    if ( h != 0 ) h = h->Clone();
    if ( h != 0 ) h->Write();
    cout << ". to get the list" << endl;
  }
  cout << endl;

  delete file;

  return 0;

}

