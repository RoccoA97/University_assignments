#include <iostream>
#include <fstream>
#include <string>

#include "EventReadFromFile.h"
#include "Event.h"

using namespace std;


// read data from file "name"
EventReadFromFile::EventReadFromFile( const string& name ) {
  file = new ifstream( name.c_str(), ios::binary );
}


EventReadFromFile::~EventReadFromFile() {
  delete file;
}


// get an event
const Event* EventReadFromFile::get() {
  return readFile();
}


// read an event
const Event* EventReadFromFile::readFile() {

  float file_x;
	float file_y;
	float file_z;
	int file_eventId;
	int file_nPart;

	Event* ev;

	if ( (*file >> file_eventId) != 0 ) {
		*file >> file_x >> file_y >> file_z;

		ev = new Event( file_eventId, file_x, file_y, file_z );
	}

	else
		return 0;


  int file_charge;
	float file_px;
	float file_py;
	float file_pz;

	*file >> file_nPart;
	for (int j = 0; j < file_nPart; ++j) {
		*file >> file_charge;
		*file >> file_px >> file_py >> file_pz;

		ev->add( file_px, file_py, file_pz, file_charge );
	}

	return ev;

}
