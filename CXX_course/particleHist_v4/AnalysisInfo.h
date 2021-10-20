#ifndef AnalysisInfo_h
#define AnalysisInfo_h

#include <string>
#include <vector>

class AnalysisInfo {

 public:

  AnalysisInfo( int argc, char* argv[] );
  ~AnalysisInfo();

  // get the list of all words
  const std::vector<std::string>& argList() const;
  // get the word coming after the word 'key'
  const std::string& value( const std::string& key ) const;
  // check if the word 'key' is present
  bool contains( const std::string& key ) const;

 private:

  // container for words
  std::vector<std::string> args;
  // value to return for keys not found
  static std::string defaultString;

  AnalysisInfo           ( const AnalysisInfo& x );
  AnalysisInfo& operator=( const AnalysisInfo& x );

};

#endif

