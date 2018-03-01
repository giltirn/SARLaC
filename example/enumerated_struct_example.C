//Demonstrate the usage of the enumerated_struct macros
#include<containers/enumerated_struct.h>

using namespace CPSfit;

//Define for a struct 'Params' of 2 double members A and B with defaults 1.0 and 2.0 respectively
#define PARAMS (Params, double, (A)(B), (1.0)(2.0))

//Define the struct outside class scope
DEF_ENUMERATED_STRUCT(PARAMS);

//Define the struct inside class scope (these won't collide)
struct Enclosed{
  DEF_ENUMERATED_STRUCT_MEMBER(Enclosed,PARAMS);
};
DEF_ENUMERATED_STRUCT_MEMBER_EXTERNAL(Enclosed,PARAMS);


int main(void){
  //Demonstrate ET, constructor, text output
  Params p(2.3,3.2);
  std::cout << "p: " << p << std::endl;

  p = 2.0 * p;
  std::cout << "2.0 * p: "<<  p << std::endl;  

  //Demonstrate parsing
  Enclosed::Params q;
  std::string qstr = "{ A=9.2 B=8.8 }";
  std::istringstream ii(qstr);
  boost::spirit::istream_iterator iii(ii);
  iii >> q;

  std::cout << "Parsed string \"" << qstr << "\" into: "<< q << std::endl; 
  
  //Demonstrate HDF5 serialization
#ifdef HAVE_HDF5
  {
    HDF5writer wr("test.hdf5");
    write(wr,q,"q");
  }
  Enclosed::Params r;
  {
    HDF5reader rd("test.hdf5");
    read(rd,r,"q");
  }
  std::cout << "Serialized to test.hdf5 and re-read: " << r << std::endl;
  assert( r == q ); //demonstrate equivalence operator
#endif

  return 0;
}
