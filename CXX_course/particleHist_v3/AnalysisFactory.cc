#include "AnalysisFactory.h"
#include "AnalysisInfo.h"

#include <iostream>

using namespace std;

AnalysisFactory::AnalysisFactory() {
}


AnalysisFactory::~AnalysisFactory() {
}


vector<AnalysisSteering*> AnalysisFactory::create( const AnalysisInfo* info ) {
  vector<AnalysisSteering*> aList;
  // loop over analysis object factories
  static map<string,AbsFactory*>* fm = factoryMap();
  map<string,AbsFactory*>::iterator iter = fm->begin();
  map<string,AbsFactory*>::iterator iend = fm->end();
  while ( iter != iend ) {
    const pair<string,AbsFactory*>& element = *iter++;
    // create analysis object if its name is listed in the command line
    if ( info->contains( element.first ) )
        aList.push_back( element.second->create( info ) );
  }
  return aList;
}


// function to add analyzer concrete analyzer factories to the factory
void AnalysisFactory::registerFactory( const string& name, AbsFactory* b ) {
  static map<string,AbsFactory*>& fm = *factoryMap();
  fm[name] = b;
  return;
}


// map to associate analyzer names with corresponding factories
std::map<std::string,AnalysisFactory::AbsFactory*>*
                     AnalysisFactory::factoryMap() {
  static map<string,AbsFactory*>* fm = new map<string,AbsFactory*>;
  return fm;
}

