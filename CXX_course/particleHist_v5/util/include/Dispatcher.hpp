//#include "Dispatcher.h"
#include "LazyObserver.h"
#include "ActiveObserver.h"

#include <iostream>

template <class T>
const T* Dispatcher<T>::last = 0;

template <class T>
Dispatcher<T>::Dispatcher() {
}

template <class T>
Dispatcher<T>::~Dispatcher() {
}

template <class T>
void Dispatcher<T>::  subscribe( ActiveObserver<T>* obs ) {
  activeObserverList()->insert( obs );
  return;
}

template <class T>
void Dispatcher<T>::unsubscribe( ActiveObserver<T>* obs ) {
  activeObserverList()->erase( obs );
  return;
}

template <class T>
void Dispatcher<T>::  subscribe( LazyObserver<T>* obs ) {
  lazyObserverList()->insert( obs );
  if ( last != 0 ) obs->lazyUpdate( *last );
  return;
}

template <class T>
void Dispatcher<T>::unsubscribe(   LazyObserver<T>* obs ) {
  lazyObserverList()->erase( obs );
  return;
}

template <class T>
void Dispatcher<T>::notify( const T& x ) {

  last = &x;

  static   std::set<  LazyObserver<T>*>* lol =   lazyObserverList();
  typename std::set<  LazyObserver<T>*>::iterator l_iter = lol->begin();
  typename std::set<  LazyObserver<T>*>::iterator l_iend = lol->end();
  while ( l_iter != l_iend ) (*l_iter++)->lazyUpdate( x );

  static   std::set<ActiveObserver<T>*>* aol = activeObserverList();
  typename std::set<ActiveObserver<T>*>::iterator a_iter = aol->begin();
  typename std::set<ActiveObserver<T>*>::iterator a_iend = aol->end();
  while ( a_iter != a_iend ) (*a_iter++)->update( x );

  return;

}

template <class T>
std::set<ActiveObserver<T>*>* Dispatcher<T>::activeObserverList() {
  static std::set<ActiveObserver<T>*>* ptr =
     new std::set<ActiveObserver<T>*>;
  return ptr;
}

template <class T>
std::set<  LazyObserver<T>*>* Dispatcher<T>::  lazyObserverList() {
  static std::set<  LazyObserver<T>*>* ptr =
     new std::set<  LazyObserver<T>*>;
  return ptr;
}

