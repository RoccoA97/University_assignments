#include <fstream>

#include "Event.h"


const Event* read( std::ifstream& iFile ) {

	float file_x;
	float file_y;
	float file_z;
	int file_eventId;
	int file_nPart;

	Event* ev;

	if ( (iFile >> file_eventId) != 0 ) {
		iFile >> file_x >> file_y >> file_z;

		ev = new Event( file_eventId, file_x, file_y, file_z );
	}

	else
		return 0;


	iFile >> file_nPart;
	for (int j = 0; j < file_nPart; ++j) {

		int file_charge;
		float file_px;
		float file_py;
		float file_pz;

		iFile >> file_charge;
		iFile >> file_px >> file_py >> file_pz;


		ev->add( file_px, file_py, file_pz, file_charge );
	}

	return ev;
}
