#ifndef Dispatcher_H
#define Dispatcher_H

#include <set>

template <class T> class ActiveObserver;
template <class T> class   LazyObserver;

template <class T>
class Dispatcher {

  friend class ActiveObserver<T>;
  friend class   LazyObserver<T>;

 public:

  static void notify( const T& x );
  static std::set<ActiveObserver<T>*>* init();

 private:

  static void   subscribe( ActiveObserver<T>* obs );
  static void unsubscribe( ActiveObserver<T>* obs );
  static void   subscribe(   LazyObserver<T>* obs );
  static void unsubscribe(   LazyObserver<T>* obs );

  static const T* last;
  static std::set<ActiveObserver<T>*>* activeObserverList();
  static std::set<  LazyObserver<T>*>*   lazyObserverList();

  Dispatcher();
  ~Dispatcher();

};

#include "Dispatcher.hpp"

#endif

