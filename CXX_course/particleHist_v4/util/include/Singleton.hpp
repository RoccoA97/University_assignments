// included by "Singleton.h"
#include <iostream>

template <class T>
bool Singleton<T>::verbose = false;

template <class T>
Singleton<T>::Singleton() {
  if ( verbose ) std::cout << "create Singleton " << this << std::endl; 
}

template <class T>
Singleton<T>::~Singleton() {
  if ( verbose ) std::cout << "delete Singleton " << this << std::endl; 
}

template <class T>
T* Singleton<T>::instance() {
  if ( verbose ) std::cout << "Singleton::instance " << std::endl; 
  // the object is created only once, the first time "instance()" is called
  static T* obj = new T;
  return obj;
}

