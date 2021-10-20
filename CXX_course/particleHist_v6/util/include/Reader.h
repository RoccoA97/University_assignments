#ifndef Reader_h
#define Reader_h

#include <iostream>
#include <fstream>
#include <string>

// interface classes to give uniform access to text and binary files
class Reader {
 public:
  // constructor and destructor
  Reader() {}
  virtual ~Reader() { delete is; }
  // input operator for generic objects
  template <class T>
  Reader& operator>>( T& x ) {
    // create wrapper and do input operation
    Wrapper<T> w( x );
    read( w );
    // return *this to allow input operators concatenation
    return *this;
  }
  // type conversion
  operator bool() { return status; }
 protected:
  // input stream and status
  std::istream* is;
  bool status;
  // visitor wrapper base for input operations
  class BaseWrapper {
   public:
    BaseWrapper() {}
    // cannot rely on overloading to select required read function,
    // concrete acceptors (i.e. reader) classes not yet defined
    virtual bool readT( Reader* r ) = 0; // read from text   file
    virtual bool readB( Reader* r ) = 0; // read from binary file
  };
 private:
  // Concrete visitor for generic input types
  template <class T>
  class Wrapper: public BaseWrapper {
   public:
    // constructor and input functions
    Wrapper( T& x ): p( &x ) {}
    virtual bool readT( Reader* r ) {
      return static_cast<bool>( *r->is >> *p );
    }
    virtual bool readB( Reader* r ) {
      return static_cast<bool>( r->is->read( reinterpret_cast<char*>( p ),
                                          sizeof( *p ) ) );
    }
   private:
    // input object pointer
    T* p;
  };
  // accept (read) function
  virtual void read( BaseWrapper& w ) = 0;
};

// concrete reader for text input files
class TxtReader: public Reader {
 public:
  TxtReader( const std::string& name = "" ) {
    if ( name != "" )
    is = new std::ifstream( name.c_str() );
    else
    is = &std::cin;
  }
 private:
  virtual void read( BaseWrapper& w ) { status = w.readT( this ); } 
};

// concrete reader for binary input files
class BinReader: public Reader {
 public:
  BinReader( const std::string& name ) {
    is = new std::ifstream( name.c_str(), std::ios::binary );
  }
 private:
  virtual void read( BaseWrapper& w ) { status = w.readB( this ); } 
};

// reader factory
class ReaderFactory {
 public:
  static Reader* create( const std::string& name,
                         std::ios::openmode mode = std::ios::in ) {
    if ( mode & std::ios::binary )
    return new BinReader( name );
    return new TxtReader( name );
  }
};

#endif

