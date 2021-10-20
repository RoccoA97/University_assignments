#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

int main( int argc, char* argv[] ) {

  string name( argv[1] );
  vector<string> type;
  vector<string> func;
  type.reserve( argc );
  func.reserve( argc );
  int iarg = 2;
  while ( iarg < argc - 1 ) {
    type.push_back( string( argv[iarg++] ) );
    func.push_back( string( argv[iarg++] ) );
  }

  int nfunc = func.size();
  int ifunc;
  for ( ifunc = 0; ifunc < nfunc; ++ ifunc ) cout << type[ifunc] << " " << func[ifunc] << endl;

  string int_name = name + ".h";
  string imp_name = name + ".cc";
  ofstream int_file( int_name.c_str() );
  ofstream imp_file( imp_name.c_str() );

  int_file << "#ifndef " << name << "_h" << endl;
  int_file << "#define " << name << "_h" << endl;
  int_file << endl;
  int_file << "class " << name << " {" << endl;
  int_file << endl;
  int_file << " public:" << endl;
  int_file << endl;
  int_file << "  "          << name << "();" << endl;
  int_file << "  virtual ~" << name << "();" << endl;
  for ( ifunc = 0; ifunc < nfunc; ++ ifunc ) {
  int_file << endl;
  int_file << "  " << type[ifunc] << " " << func[ifunc] << "();" << endl;
  }
  int_file << endl;
  int_file << " private:" << endl;
  int_file << endl;
  int_file << "  " << name << "           ( const " << name << "& x );" << endl;
  int_file << "  " << name << "& operator=( const " << name << "& x );" << endl;
  int_file << endl;
  int_file << "};" << endl;
  int_file << endl;
  int_file << "#endif" << endl;
  int_file << endl;

  imp_file << "#include \"" << name << ".h\"" << endl;
  imp_file << endl;
  imp_file << endl;
  imp_file << name << "::" << name << "() {" << endl;
  imp_file << "}" << endl;
  imp_file << endl;
  imp_file << endl;
  imp_file << name << "::~" << name << "() {" << endl;
  imp_file << "}" << endl;
  for ( ifunc = 0; ifunc < nfunc; ++ ifunc ) {
  imp_file << endl;
  imp_file << endl;
  imp_file << type[ifunc] << " " << name << "::" << func[ifunc] << "() {"
           << endl;
  imp_file << "}" << endl;
  }
  imp_file << endl;

  return 0;

}

