#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "Event.h"
#include "EventSource.h"
#include "SourceFactory.h"
#include "AnalysisInfo.h"
#include "AnalysisSteering.h"
#include "AnalysisFactory.h"

using namespace std;


int main( int argc, char* argv[] ) {

  // store command line parameters
  AnalysisInfo* info = new AnalysisInfo( argc, argv );

  // create data source
  EventSource* es = SourceFactory::create( info );

  // create a list of analyzers
  vector<AnalysisSteering*> aList = AnalysisFactory::create( info );

  // variables to loop over analyzers
  int l = aList.size();
  int i;

  // initialize all analyzers
  for ( i = 0; i < l; ++ i ) aList[i]->beginJob();

  // loop over events
  const Event* ev;
  while ( ( ev = es->get() ) != 0 ) {
    for ( i = 0; i < l; ++ i ) aList[i]->process( *ev );
    delete ev;
  }

  // finalize all analyzers
  for ( i = 0; i < l; ++ i ) {
    aList[i]->endJob();
    cout << endl;
  }

  return 0;

}
