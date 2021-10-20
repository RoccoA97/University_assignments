#include <iostream>

#include "Event.h"


// Dump data to the terminal
void dump( const Event& ev ) {

	std::cout << ev.eventNumber() << std::endl
						<< ev.getX()				<< ' '
						<< ev.getY()				<< ' '
						<< ev.getZ()				<< std::endl
						<< ev.nParticles() 	<< std::endl;

	// Loop over ev.nParticles() particles
	for ( unsigned int i = 0; i < ev.nParticles(); ++i ) {
		std::cout << (ev.particle(i))->charge << ' '
							<< (ev.particle(i))->px 		<< ' '
							<< (ev.particle(i))->py 		<< ' '
							<< (ev.particle(i))->pz 		<< std::endl;
	}


	return;
}
