//#include "ActiveObserver.h"
#include "Dispatcher.h"
#include <iostream>

template <class T>
ActiveObserver<T>::ActiveObserver() {
  Dispatcher<T>::subscribe( this );
}

template <class T>
ActiveObserver<T>::~ActiveObserver() {
  Dispatcher<T>::unsubscribe( this );
}

