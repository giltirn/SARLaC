#include<sstream>
#include<array>
#include<parser/parser.h>

using namespace SARLaC;

#define TEST_MEMBERS \
  ( int, a) \
  ( std::string, b)

struct Test{
  GENERATE_MEMBERS(TEST_MEMBERS);
};
GENERATE_PARSER(Test, TEST_MEMBERS);

#define TEST2_MEMBERS \
  ( Test, a)

struct Test2{
  GENERATE_MEMBERS(TEST2_MEMBERS);
};
GENERATE_PARSER(Test2, TEST2_MEMBERS);

GENERATE_ENUM_AND_PARSER(EnumTest, (A)(B) );

int main(const int argc, const char* argv[]){


  //Test array parse
  {
    std::stringstream is;
    is << "(1,2,3)";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    std::array<int,3> v;
    f >> v;

    assert(v[0] == 1 && v[1] == 2 && v[2] == 3);
  }

  //Test array parse error handling
  {
    std::stringstream is;
    is << "(1,2)";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    std::array<int,3> v;
    bool failed = false;
    try{
      std::cout << "This should fail:\n";
      f >> v;
    }catch(...){
      std::cout << "Got expected failure!\n";
      failed = true;
    }
    assert(failed);
  }
    
  //Test struct parse
  {
    std::stringstream is;
    is << "{ a=5 b=\"hello\" }";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    Test v;
    f >> v;
    std::cout << v.a << " " << v.b << std::endl;
    assert(v.a == 5 && v.b == "hello");
  }

  //Test struct error handling
  {
    std::stringstream is;
    is << "{ a=5.4 b=\"hello\" }";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    Test v;
    bool failed = false;
    try{
      std::cout << "This should fail:\n";
      f >> v;
    }catch(std::exception const &ex){
      std::cerr << ex.what() << std::endl;
      std::cout << "Got expected failure!\n";
      failed = true;
    }
    assert(failed);
  }


  //Test compound struct parse
  {
    std::stringstream is;
    is << "{ a = { a=5 b=\"hello\" } }";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    Test2 v;
    f >> v;
    std::cout << v.a.a << " " << v.a.b << std::endl;
    assert(v.a.a == 5 && v.a.b == "hello");
  }

  //Test struct error handling
  {
    std::stringstream is;
    is << "{ a = { q=5 b=\"hello\" } }";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    Test2 v;
    bool failed = false;
    try{
      std::cout << "This should fail:\n";
      f >> v;
    }catch(std::exception const &ex){
      std::cerr << ex.what() << std::endl;
      std::cout << "Got expected failure!\n";
      failed = true;
    }
    assert(failed);
  }

  //Test enum parser
  {
    std::stringstream is;
    is << "A";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    EnumTest v;
    f >> v;
    std::cout << v << std::endl;
    assert(v == EnumTest::A);
  }
  
  //Test enum error handling
  {
    std::stringstream is;
    is << "";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    EnumTest v;
    bool failed = false;
    try{
      std::cout << "This should fail:\n";
      f >> v;
    }catch(std::exception const &ex){
      std::cerr << ex.what() << std::endl;
      std::cout << "Got expected failure!\n";
      failed = true;
    }
    assert(failed);
  }

  //Test std::pair parse
  {
    std::stringstream is;
    is << "{ A, 2.3 }";
    is >> std::noskipws;
    boost::spirit::istream_iterator f(is);

    std::pair<char,double> v;
    f >> v;
    std::cout << v.first << " " << v.second << std::endl;
    assert(v.first == 'A' && v.second == 2.3);
  }

  //Test std::map rule
  {
    std::stringstream ss; 
    ss << "{ \"Hello\":1,  \"World\":2 }";
    ss >> std::noskipws;
    boost::spirit::istream_iterator f(ss);
    
    std::map<std::string,int> v;
    f >> v;
    for(auto it=v.begin(); it !=v.end(); it++){
      std::cout << it->first << " : " << it->second << std::endl;
    }
    auto it = v.find("Hello");
    assert(it != v.end() && it->second == 1);
    it = v.find("World");
    assert(it != v.end() && it->second == 2);
  } 
    

  

  std::cout << "Test complete with no (unexpected!) errors" << std::endl;
  return 0;
}
