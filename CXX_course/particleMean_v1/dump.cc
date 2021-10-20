#include <iostream>

#include "Event.h"


// Dump data to the terminal
void dump( const Event& ev ) {

	std::cout << ev.eventId << std::endl
						<< ev.x				<< ' '
						<< ev.y				<< ' '
						<< ev.z				<< std::endl
						<< ev.nPart 	<< std::endl;

	// Loop over ev.nPart particles
	for ( int i = 0; i < ev.nPart; ++i ) {
		std::cout << ((ev.pt)[i])->charge << ' '
							<< ((ev.pt)[i])->px 		<< ' '
							<< ((ev.pt)[i])->py 		<< ' '
							<< ((ev.pt)[i])->pz 		<< std::endl;
	}


	return;
}
