#ifndef ActiveObserver_H
#define ActiveObserver_H

template <class T>
class ActiveObserver {

 public:

  ActiveObserver();
  virtual ~ActiveObserver();

  virtual void update( const T& x ) = 0;

 private:

  ActiveObserver           ( const ActiveObserver& x );
  ActiveObserver& operator=( const ActiveObserver& x );

};

#include "ActiveObserver.hpp"

#endif

