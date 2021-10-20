#ifndef EventSource_h
#define EventSource_h

class Event;

class EventSource {

 public:

  EventSource();
  virtual ~EventSource();

  // get events
  virtual void run();


 private:

  EventSource           ( const EventSource& x );
  EventSource& operator=( const EventSource& x );

  // get an event
  virtual const Event* get() = 0;

};

#endif
