#include<sstream>
#include<array>
#include<parser/parser.h>

using namespace CPSfit;


int main(const int argc, const char* argv[]){
  //Test array parse
  std::stringstream is;
  is << "(1,2,3)";

  is >> std::noskipws;
  boost::spirit::istream_iterator f(is);

  std::array<int,3> v;
  f >> v;

  for(int i=0;i<3;i++) std::cout << v[i] << " ";
  std::cout << std::endl;

  assert(v[0] == 1 && v[1] == 2 && v[2] == 3);
  
  return 0;
}
