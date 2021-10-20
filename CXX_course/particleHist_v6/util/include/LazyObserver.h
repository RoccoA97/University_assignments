#ifndef LazyObserver_H
#define LazyObserver_H

template <class T>
class LazyObserver {

 public:

  LazyObserver();
  virtual ~LazyObserver();

  virtual void lazyUpdate( const T& x );

 protected:

  virtual void update( const T& x ) = 0;
  virtual void check();

 private:

  bool upToDate;
  bool updating;
  const T* last;

  LazyObserver           ( const LazyObserver& x );
  LazyObserver& operator=( const LazyObserver& x );

};

#include "LazyObserver.hpp"

#endif

