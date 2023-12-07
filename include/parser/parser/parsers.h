#ifndef _SARLAC_PARSERS_H___
#define _SARLAC_PARSERS_H___

#include <vector>
#include <map>
#include <boost/spirit/home/x3/support/utility/error_reporting.hpp>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils/type_info.h>
#include<utils/template_wizardry/text_io.h>
#include<parser/parser/parser_tools.h>

SARLAC_START_NAMESPACE

//Bindings of types to parsers. Create overloads automatically or by hand
namespace parsers{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;

  struct error_handler
  {
    template <typename Iterator, typename Exception, typename Context>
    x3::error_handler_result on_error(
				      Iterator& first, Iterator const& last
				      , Exception const& x, Context const& context)
    {
      auto& error_handler = x3::get<x3::error_handler_tag>(context).get();
      std::string message = "Error! Expecting: " + x.which() + " here:";
      error_handler(x.where(), message);
      //return x3::error_handler_result::fail;
      return x3::error_handler_result::rethrow;
    }
  };


  template<typename T>
  struct parser{};

  DEF_X3_PARSER(char);
  DEF_X3_PARSER(int);
  DEF_X3_PARSER(double);
  DEF_X3_PARSER(bool);

#define STRIP_PARENS(...) __VA_ARGS__

template<typename T>
struct getAttribute{
  typedef decltype(std::declval<T>().parse) rule;
  typedef typename rule::attribute_type attr;
  typedef typename rule::id id;
};

#ifdef OLD_BOOST_X3
#define DEF_BASIC_PARSE_RULE_PEXPL(templ, type)				\
  template <STRIP_PARENS templ typename Iterator, typename Context, typename Attribute> \
  bool parse_rule(decltype(std::declval<parser<STRIP_PARENS type> >().parse) rule_ , Iterator& first, Iterator const& last , Context const& context, Attribute& attr){ \
    using boost::spirit::x3::unused;					\
    static parser<STRIP_PARENS type> inst;				\
    static auto const parse_def = inst.get_def();			\
    static auto const def_ = (inst.parse = parse_def);			\
    return def_.parse(first, last, context, unused, attr);		\
  }
#else
#define DEF_BASIC_PARSE_RULE_PEXPL(templ, vtype, parsertype)				\
  template <STRIP_PARENS templ typename Iterator, typename Context>	\
  inline bool parse_rule( ::boost::spirit::x3::detail::rule_id<typename decltype(std::declval<STRIP_PARENS parsertype>().parse)::id>, Iterator& first, Iterator const& last , Context const& context, STRIP_PARENS vtype & attr) { \
    static STRIP_PARENS parsertype inst;					\
    static auto const parse_def = inst.get_def();				\
    static auto const def_ = (inst.parse = parse_def);			\
    typedef decltype(inst.parse) rule_t;				\
									\
    return ::boost::spirit::x3::detail ::rule_parser<typename rule_t::attribute_type, typename rule_t::id, true> ::call_rule_definition( def_, inst.parse.name , first, last, context, attr , ::boost::mpl::bool_<rule_t::force_attribute>()); \
  }

#endif
#define DEF_BASIC_PARSE_RULE(templ, vtype) DEF_BASIC_PARSE_RULE_PEXPL(templ, vtype, (parser<STRIP_PARENS vtype>))

  //Vector parser
  //For generic vector parse we can hack the parser interface to allow templating the underlying type by static instantiating the rule definition inside the parse_rule itself.
  template<typename U>
  class vector_T_rule: error_handler{};

  template<typename T>
  struct parser< std::vector<T> >{
    typedef std::vector<T> type;
    x3::rule<class vector_T_rule<T>, std::vector<T> > const parse;
    parsers::parser<T> elem_parser;
    
    parser(): parse( typeid(std::vector<T>).name() ){}

    inline decltype(auto) get_def() const{
      return x3::char_('(') >> *( elem_parser.parse[parser_tools::push_back] >> *(',' >> elem_parser.parse[parser_tools::push_back]) ) > ')';
    }
  };
  DEF_BASIC_PARSE_RULE((typename T,), (std::vector<T>));

  //Array parser
  template<typename T, std::size_t S>
  class array_T_rule: error_handler{};
  
  template<typename T, std::size_t I, std::size_t S>
  struct gen_array_rule_recurse{
    template<typename P>
    static inline decltype(auto) recurse(const parsers::parser<T> &elem_parser, const P &parser){ 
      auto NP = parser > elem_parser.parse[parser_tools::set_elem<T>(I-1)] > ','; 
      return gen_array_rule_recurse<T,I+1,S>::recurse(elem_parser,NP);
    }
  };
  template<typename T, std::size_t S>
  struct gen_array_rule_recurse<T,S,S>{
    template<typename P>
    static inline decltype(auto) recurse(const parsers::parser<T> &elem_parser, const P &parser){ 
      return parser > elem_parser.parse[parser_tools::set_elem<T>(S-1)]; 
    }
  };

  template<typename T, std::size_t S>
  struct parser< std::array<T,S> >{
    x3::rule<class array_T_rule<T,S>, std::array<T,S> > const parse;
    parsers::parser<T> elem_parser;
    
    parser(): parse( typeid(std::array<T,S>).name() ){}

    inline decltype(auto) get_def() const{
      auto NP = x3::char_('(');
      auto NP2 = gen_array_rule_recurse<T,1,S>::recurse(elem_parser,NP);
      return NP2 > ')';
    }
  };
  DEF_BASIC_PARSE_RULE((typename T, size_t S,), (std::array<T,S>));

  //Pair parser
  template<typename T, typename U>
  class pair_T_U_rule: error_handler{};

  template<typename T, typename U>
  struct parser< std::pair<T,U> >{
    x3::rule<class pair_T_U_rule<T,U>, std::pair<T,U> > const parse;
    parsers::parser<T> T_parser;
    parsers::parser<U> U_parser;
    
    parser(): parse( typeid(std::pair<T,U>).name() ){}

    inline decltype(auto) get_def() const{
      return x3::char_('{') > T_parser.parse[parser_tools::member_set_equals<std::pair<T,U>,T,&std::pair<T,U>::first>()] > ','
	> U_parser.parse[parser_tools::member_set_equals<std::pair<T,U>,U,&std::pair<T,U>::second>()] > '}';
    }
  };
  DEF_BASIC_PARSE_RULE((typename T, typename U,), (std::pair<T,U>));

  
  //String parser
  class string_rule_: error_handler{};
  auto const quoted_string = x3::lexeme['"' > *(x3::char_ - '"') > '"'];
  x3::rule<string_rule_, std::string> const string_rule = "string_rule";
  auto const string_rule_def = quoted_string[parser_tools::set_equals];
  BOOST_SPIRIT_DEFINE(string_rule);
  DEF_CUSTOM_PARSER(std::string, string_rule);

  //Complex parser
  class complexD_: error_handler{};
  auto const complexD_rule_def = x3::no_skip[*x3::lit(' ') > x3::double_[parser_tools::set_reim<double,0>()] > x3::lit(' ') > x3::double_[parser_tools::set_reim<double,1>()]];
  x3::rule<complexD_, std::complex<double> > const complexD_rule = "complexD_rule";
  BOOST_SPIRIT_DEFINE(complexD_rule);
  DEF_CUSTOM_PARSER(std::complex<double>, complexD_rule);

  //std::map parser
  template<typename T, typename U>
  class map_pair_T_U_rule: error_handler{};

  template<typename T, typename U>
  struct map_pair_T_U_parser{
    x3::rule<class map_pair_T_U_rule<T,U>, std::pair<T,U> > const parse;
    parsers::parser<T> T_parser;
    parsers::parser<U> U_parser;
    
    map_pair_T_U_parser(): parse( typeid(std::pair<T,U>).name() ){}

    inline decltype(auto) get_def() const{
      return T_parser.parse[parser_tools::member_set_equals<std::pair<T,U>,T,&std::pair<T,U>::first>()] > ':'
  	> U_parser.parse[parser_tools::member_set_equals<std::pair<T,U>,U,&std::pair<T,U>::second>()];
    }
  };
  DEF_BASIC_PARSE_RULE_PEXPL((typename T, typename U,), (std::pair<T,U>), (map_pair_T_U_parser<T,U>));

  template<typename T, typename U>
  class map_T_U_rule: error_handler{};

  template<typename T, typename U>
  struct parser< std::map<T,U> >{
    x3::rule<class map_T_U_rule<T,U>, std::map<T,U> > const parse;
    map_pair_T_U_parser<T,U> elem_parser;

    parser(): parse( typeid(std::map<T,U>).name() ){}

    struct map_insert_elem{
      template <typename Context>
      void operator()(Context const& ctx) const{
	reinterpret_cast<std::map<T,U> &>(_val(ctx)).insert(_attr(ctx));
      }
    };

    inline decltype(auto) get_def() const{
      return x3::char_('{') >> *( elem_parser.parse[map_insert_elem()] >> *(',' >> elem_parser.parse[map_insert_elem()]) ) > '}';
    }
  };
  DEF_BASIC_PARSE_RULE((typename T, typename U,), (std::map<T,U>));
  
  template<typename T>
  struct parser_instance{
    inline static const parser<T> &get(){ 
      static parser<T> p;
      return p;
    }
  };
    

};


//Make sure the underlying stream has std::noskipws enabled or it won't parse correctly
template<typename T, typename std::enable_if<hasParseMember<parsers::parser<T> >::value, int>::type = 0>
boost::spirit::istream_iterator & operator>>(boost::spirit::istream_iterator &f, T &s){ \
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;

  boost::spirit::istream_iterator l;
  parsers::parser<T> vp;

  using boost::spirit::x3::with;
  using boost::spirit::x3::error_handler_tag;
  using error_handler_type = boost::spirit::x3::error_handler<boost::spirit::istream_iterator>;

  error_handler_type error_handler(f, l, std::cerr);

  auto const parser =
    with<error_handler_tag>(std::ref(error_handler))
    [
     vp.parse
     ];

  bool r = x3::phrase_parse(f, l, parser, ascii::space, s);

  if(!r){
    throw std::runtime_error(std::string("Parsing of type ") + printType<T>() +std::string(" failed\n"));
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





SARLAC_END_NAMESPACE
#endif
