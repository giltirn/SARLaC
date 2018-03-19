#ifndef _PARSER_H_
#define _PARSER_H_

#include<config.h>
#include<utils/macros.h>
#include<parser/parser/parsers.h>

CPSFIT_START_NAMESPACE

//Combine all the above the specify a parser *outside the CPSfit namespace*
//Usage: GENERATE_PARSER( MY_STRUCT_NAME ,  (MY_TYPE1, MY_MEMBER1)(MY_TYPE2, MY_MEMBER2).... )
#define GENERATE_PARSER(structname, structmembers) _GENERATE_PARSER(structname, structmembers)

//Same as the above but with the grammar namespace manually specified
#define GENERATE_PARSER_GM(structname, grammar, structmembers) _GENERATE_PARSER_GM(structname, grammar, structmembers)

//Convenience call for generating the class members using the same interface as for the parser (optional)
//Usage: GENERATE_MEMBERS( (MY_TYPE1, MY_MEMBER1)(MY_TYPE2, MY_MEMBER2).... )
#define GENERATE_MEMBERS(structmembers) _GENERATE_MEMBERS(structmembers) 

//Define a parser for an enum
//Usage: GENERATE_ENUM_PARSER( MY_ENUM NAME, (MY_ELEM1)(MY_ELEM2)(MY_ELEM3)... )
#define GENERATE_ENUM_PARSER(enumname, enummembers) _GENERATE_ENUM_PARSER(enumname, enummembers) 

//Define an enum and a parser to go along with it
//Usage: GENERATE_ENUM_AND_PARSER( MY_ENUM NAME, (MY_ELEM1)(MY_ELEM2)(MY_ELEM3)......   )
#define GENERATE_ENUM_AND_PARSER(enumname, enummembers) _GENERATE_ENUM_AND_PARSER(enumname, enummembers) 
  
//Parse a file to a struct. Struct type must have a parser created using the above macros
template<typename T, typename std::enable_if<hasParseMember<parsers::parser<T> >::value, int>::type = 0 >
void parse(T &s, const std::string &filename){
  std::ifstream is(filename.c_str());
  if(is.fail()) error_exit(std::cout << "parse(T&,const std::string&) error when opening file " << filename << ": " << strerror(errno) << std::endl);
  is >> std::noskipws;
  boost::spirit::istream_iterator f(is);
  f >> s;
}

CPSFIT_END_NAMESPACE

#endif
