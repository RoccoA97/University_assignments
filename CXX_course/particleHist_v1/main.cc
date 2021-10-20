#include <vector>
#include <string>
#include <sstream>

#include "Event.h"
#include "EventReadFromFile.h"
#include "EventSim.h"
#include "AnalysisSteering.h"
#include "ParticleMass.h"
#include "EventDump.h"

using namespace std;


int main( int argc, char* argv[] ) {

  // create data source
  EventSource* es;
  const string type = argv[1];
  if ( type == "input" ) {
    const string name = argv[2];
    es = new EventReadFromFile( name );
  }

  else if ( type == "sim" ) {
    const string nevt = argv[2];
    const string seed = ( argc > 3 ? argv[3] : "1" );
    stringstream sstr;
    unsigned int n;
    sstr.str( nevt );
    sstr >> n;
    sstr.clear();
    unsigned int s;
    sstr.str( seed );
    sstr >> s;
    es = new EventSim( n, s );
  }

  else {
    cout << "invalid keyword" << endl;
    return 0;
  }


  // create a list of analyzers
  vector<AnalysisSteering*> aList;

  // create object to dump event
  // and store into list of analyzers
  aList.push_back( new EventDump );

  // create object to compute mean and rms masses
  // and store into list of analyzers
  aList.push_back( new ParticleMass );

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
