#ifndef TFileProxy_h
#define TFileProxy_h

// Proxy to handle repeated open and close of the same TFile
//
// If a TFile is opened several times with the same mode, at the first
// access it's opened, or the operation is refused (for existing files
// opened with "NEW" or "CREATE" mode), but at the following accesses 
// an undesidered behaviour can appear:
// "NEW" or "CREATE": after closing the file cannot be opened anymore
//                    with the same mode
// "RECREATE"       : after closing the file can be opened again, but
//                    it's overwritten
// On the other side "UPDATE" mode allows to open the file again without 
// overwriting, but that prevents the checks on the file existence
// made with "CREATE" and "RECREATE" modes.
//
// A proxy can solve this problem, taking care of using the specified
// access mode the first time the file is opened, and using "UPDATE"
// mode at the following accesses (but for readonly files)


#include "TFile.h"

#include <string>
#include <map>

class TFileProxy {

 public:

  // constructor
  TFileProxy( const std::string& name, const std::string& mode ):
   fName( name ) {
    // look for file in map
    static std::map<std::string,TFile*>& m = fMap();
    std::map<std::string,TFile*>::iterator iter = m.find( name );
    // if not found open with specified mode
    if ( iter == m.end() ) {
      if ( !( m[name] = new TFile( name.c_str(), mode.c_str() ) )->IsOpen() )
              m.erase( name );
    }
    // if found reopen the file
    else {
      // reopen file for writing
      if ( ( mode != "" ) && ( mode != "READ" ) )
           iter->second = new TFile( name.c_str(), "UPDATE" );
      // reopen readonly file
      else iter->second = new TFile( name.c_str() );
      // set current directory to file
      iter->second->cd();
    }
  }

  // destructor
  ~TFileProxy() {
    Close();
  }

  // forward functions

  bool cd() {
    static std::map<std::string,TFile*>& m = fMap();
    return m[fName]->cd();
  }

  void Close() {
    static TFile* fnp = 0;
    static std::map<std::string,TFile*>& m = fMap();
    std::map<std::string,TFile*>::iterator iter = m.find( fName );
    TFile*& file = ( iter == m.end() ? fnp : iter->second );
    if ( file == 0 ) return;
    file->Close();
    delete file;
    file = 0;
    return;
  }

 private:

  // file name
  std::string fName;

  // map of names to TFile
  static std::map<std::string,TFile*>& fMap() {
    static std::map<std::string,TFile*> m;
    return m;
  }

  // dummy copy constructor and assignment to prevent unadvertent copy
  TFileProxy           ( const TFileProxy& x );
  TFileProxy& operator=( const TFileProxy& x );

};

#endif

