#ifndef TDirectory_h
#define TDirectory_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

class TH1F;

class TDirectory {

  friend class TH1F;

 public:

  TDirectory();
  virtual ~TDirectory();

  void cd();
  std::ios* file() const;

  TH1F* Get( const char* name );
  std::vector<std::string> histoList();

 protected:

  std::ios* fp;
  std::map<std::string,TH1F*> hmap;

 private:

  virtual void update();

  TDirectory           ( const TDirectory& x );
  TDirectory& operator=( const TDirectory& x );

};

extern TDirectory* gDirectory;

#endif

