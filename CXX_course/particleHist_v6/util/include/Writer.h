#ifndef Writer_h
#define Writer_h

#include <iostream>
#include <fstream>
#include <string>

// interface classes to give uniform access to text and binary files
class Writer {
 public:
  // constructor and destructor
  Writer() {}
  virtual ~Writer() {
    if ( dynamic_cast<std::ofstream*>( os ) != 0 ) delete os;
  }
  // output operator for ostream functions (e.g. "endl")
  Writer& operator<<( std::ostream& ( *pf )( std::ostream& ) ) {
    *os << pf;
    return *this;
  }
  // output operator for generic objects
  template <class T>
  Writer& operator<<( const T& t ) {
    Wrapper<T> w( t );
    write( w );
    return *this;
  }
  // cleaner class to filter output instructions not suited
  // for binary files
  class BaseCleaner {
   public:
    virtual ~BaseCleaner() {}
    virtual bool write( Writer& w ) = 0;
  };
  class Cleaner {
   public:
    typedef std::ostream& ( *F )( std::ostream& );
    template <class T>
    Cleaner( const T& t ): p( new FormatCleaner<T>(  t ) ) {}
    Cleaner( F pf       ): p( new CallCleaner     ( pf ) ) {}
    ~Cleaner() { delete p; }
    BaseCleaner* get() const { return p; }
   private:
    Cleaner             ( const Cleaner& c );
    Cleaner& operator = ( const Cleaner& c );
    BaseCleaner* p;
  };
  Writer& operator<<( const Cleaner& c ) {
    write( *c.get() );
    return *this;
  }
  // type conversion
  operator bool() { return status; }
 protected:
  // output stream
  std::ostream* os;
  bool status;
  // visitor wrapper base for output operations
  class BaseWrapper {
   public:
    BaseWrapper() {}
    // cannot rely on overloading to select required write function,
    // concrete acceptors (i.e. writer) classes not yet defined
    virtual bool writeT( std::ostream* os ) = 0;
    virtual bool writeB( std::ostream* os ) = 0;
  };
 private:
  // Concrete visitor for generic output types
  template <class T>
  class Wrapper: public BaseWrapper {
   public:
    // constructor and output functions
    Wrapper( const T& t ): p( &t ) {}
    virtual bool writeT( std::ostream* os ) {
      return static_cast<bool>( *os << *p );
    }
    virtual bool writeB( std::ostream* os ) {
      return static_cast<bool>(  os->write( reinterpret_cast<const char*>( p ),
                                            sizeof( *p ) ) );
    }
   private:
    // output object pointer
    const T* p;
  };
  // Concrete cleaner for generic input types
  template<class T>
  class FormatCleaner: public BaseCleaner {
   public:
    FormatCleaner( const T& t ): p( &t ) {}
    virtual bool write( Writer& w ) { return w << *p; }
  private:
    const T* p;
  };
  // Concrete cleaner for ostream functions
  class CallCleaner: public BaseCleaner {
   public:
   CallCleaner( Cleaner::F x ):f( x ) {}
    virtual bool write( Writer& w ) { return w << f; }
  private:
    Cleaner::F f;
  };
  // accept (write) function
  virtual void write( BaseWrapper& w ) = 0;
  virtual void write( BaseCleaner& c ) = 0;
};

// concrete writer for text output files
class TxtWriter: public Writer {
 public:
  TxtWriter( const std::string& name = "" ) {
    if ( name != "" )
    os = new std::ofstream( name.c_str() );
    else
    os = &std::cout;
  }
 private:
  virtual void write( BaseWrapper& w ) { status = w.writeT( os ); } 
  virtual void write( BaseCleaner& c ) { status = c.write( *this ); } 
};

// concrete writer for binary output files
class BinWriter: public Writer {
 public:
  BinWriter( const std::string& name ) {
    os = new std::ofstream( name.c_str(), std::ios::binary );
  }
 private:
  virtual void write( BaseWrapper& w ) { status = w.writeB( os ); } 
  virtual void write( BaseCleaner& c ) {} 
};

// writer factory
class WriterFactory {
 public:
  static Writer* create( const std::string& name,
                         std::ios::openmode mode = std::ios::out ) {
    if ( mode & std::ios::binary )
    return new BinWriter( name );
    return new TxtWriter( name );
  }
};

#endif

