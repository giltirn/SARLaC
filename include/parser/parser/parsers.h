#ifndef _CPSFIT_PARSERS_H___
#define _CPSFIT_PARSERS_H___

#include<config.h>
#include<utils/macros.h>
#include<parser/parser/parser_tools.h>

CPSFIT_START_NAMESPACE

//Bindings of types to parsers. Create overloads automatically or by hand
namespace parsers{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;

  template<typename T>
  struct parser{};

  DEF_X3_PARSER(int);
  DEF_X3_PARSER(double);
  DEF_X3_PARSER(bool);

  //For generic vector parse we can hack the parser interface to allow templating the underlying type by static instantiating the rule definition inside the parse_rule itself.
  template<typename U>
  struct vector_T_rule{};

  template<typename T>
  struct parser< std::vector<T> >{
    x3::rule<class vector_T_rule<T>, std::vector<T> > const parse;
    parsers::parser<T> elem_parser;
    
    parser(): parse( typeid(std::vector<T>).name() ){}

    inline decltype(auto) get_def() const{
      return x3::char_('(') >> *( elem_parser.parse[parser_tools::push_back] >> *(',' >> elem_parser.parse[parser_tools::push_back]) ) > ')';
    }
  };
  template <typename T, typename Iterator, typename Context, typename Attribute>
  bool parse_rule(x3::rule<class vector_T_rule<T>, std::vector<T> > const rule_ , Iterator& first, Iterator const& last , Context const& context, Attribute& attr){
    using boost::spirit::x3::unused;
    static parser<std::vector<T> > inst;
    static auto const parse_def = inst.get_def();
    static auto const def_ = (inst.parse = parse_def);
    return def_.parse(first, last, context, unused, attr);
  }
  
  auto const quoted_string = x3::lexeme['"' > *(x3::char_ - '"') > '"'];
  x3::rule<struct string_rule_, std::string> const string_rule = "string_rule";
  auto const string_rule_def = quoted_string[parser_tools::set_equals];
  BOOST_SPIRIT_DEFINE(string_rule);
  DEF_CUSTOM_PARSER(std::string, string_rule);

  auto const complexD_rule_def = x3::no_skip[*x3::lit(' ') > x3::double_[parser_tools::set_reim<double,0>()] > x3::lit(' ') > x3::double_[parser_tools::set_reim<double,1>()]];
  x3::rule<struct complexD_, std::complex<double> > const complexD_rule = "complexD_rule";
  BOOST_SPIRIT_DEFINE(complexD_rule);
  DEF_CUSTOM_PARSER(std::complex<double>, complexD_rule);
};


//Make sure the underlying stream has std::noskipws enabled or it won't parse correctly
template<typename T, typename std::enable_if<hasParseMember<parsers::parser<T> >::value, int>::type = 0>
boost::spirit::istream_iterator & operator>>(boost::spirit::istream_iterator &f, T &s){ \
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;

  boost::spirit::istream_iterator l;
  parsers::parser<T> vp;

  bool r = x3::phrase_parse(f, l, vp.parse, ascii::space, s);

  if(!r){
    throw std::runtime_error("Parsing of type " BOOST_PP_STRINGIZE(NAME) " failed\n");
  }
  return f;
}

//Example hand-written parser for simple struct
/*
struct V{
  std::vector<int> v;
};
std::ostream & operator<<(std::ostream &os, const V &v){
  for(int i=0;i<v.v.size();i++) os << v.v[i] << std::endl;
  return os;
}

namespace V_parser{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;

  parsers::parser<std::vector<int> > vector_int_;
  
  x3::rule<class V_v, std::vector<int> > const V_v = "V_v";
  auto const V_v_def = x3::char_('v') >> '=' >> vector_int_.parse[parser_tools::set_equals];

  x3::rule<class V_, V> const V_ = "V_";
  auto const V__def = V_v[parser_tools::member_set_equals<V,std::vector<int>,&V::v>()];
  
  BOOST_SPIRIT_DEFINE(V_, V_v);
};

namespace parsers{
  template<>
  struct parser<V>{
    decltype( V_parser::V_ ) &parse;
    parser(): parse(V_parser::V_){}
  };

};

int main(void){
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  {
    std::string str = "v = (9,8,7)";
    V v;
    parsers::parser<V> vp;
    
    bool r = x3::phrase_parse(str.begin(), str.end(), vp.parse, ascii::space, v);  //V_parser::V_

    std::cout << "V success ? " << r << " val: " << v << std::endl;
  }
  return 0;
}
*/





CPSFIT_END_NAMESPACE
#endif