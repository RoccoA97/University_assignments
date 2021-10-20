#include "EventDump.h"
#include "Event.h"
#include "AnalysisFactory.h"

#include "util/include/ActiveObserver.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>


// concrete factory to create an EventDump analyzer
class EventDumpFactory: public AnalysisFactory::AbsFactory {
 public:
  // assign "dump" as name for this analyzer and factory
  EventDumpFactory(): AnalysisFactory::AbsFactory( "dump" ) {}
  // create an EventDump when builder is run
  virtual AnalysisSteering* create( const AnalysisInfo* info ) {
    return new EventDump( info );
  }
};
// create a global EventDumpFactory, so that it is created and registered
// before main execution start:
// when the AnalysisFactory::create function is run,
// an EventDumpFactory will be available with name "dump".
static EventDumpFactory ed;


EventDump::EventDump( const AnalysisInfo* info ):
  AnalysisSteering( info ) {
}


EventDump::~EventDump() {
}


// function to be called at execution start
void EventDump::beginJob() {
  return;
}


// function to be called at execution end
void EventDump::endJob() {
  return;
}


// function to be called for each event
void EventDump::update( const Event& ev ) {


  std::cout << ev.eventNumber() << std::endl
            << ev.getX()        << ' '
            << ev.getY()        << ' '
            << ev.getZ()        << std::endl
            << ev.nParticles()  << std::endl;

	for (unsigned int i = 0; i < (ev.nParticles()); ++i) {

		std::cout << (ev.particle(i))->charge << ' '
              << (ev.particle(i))->px     << ' '
              << (ev.particle(i))->py     << ' '
              << (ev.particle(i))->pz     << std::endl;
	}


  return;

}
