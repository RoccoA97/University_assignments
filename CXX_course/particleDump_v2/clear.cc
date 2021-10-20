#include "Event.h"


// Delete dynamic objects
// when they are not needed anymore
void clear( const Event* ev ) {

	delete[] 	ev->pt;
	delete 		ev;

	return;
}
