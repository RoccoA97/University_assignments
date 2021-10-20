#ifndef EventReadFromFile_h
#define EventReadFromFile_h

#include <iostream>
#include <fstream>
#include <string>

#include "EventSource.h"

class Event;


class EventReadFromFile: public EventSource {

 public:

  // read data from file "name"
  EventReadFromFile( const std::string& name );
  virtual ~EventReadFromFile();

  // get an event
  virtual const Event* get();

 private:

  std::ifstream* file;

  // read and event
  const Event* readFile();

  EventReadFromFile           ( const EventReadFromFile& x );
  EventReadFromFile& operator=( const EventReadFromFile& x );

};

#endif
