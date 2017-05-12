#ifndef _PARSER_H_
#define _PARSER_H_
//#define BOOST_SPIRIT_X3_DEBUG
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <utils.h>

namespace parser_tools{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  
  //Use these as lambdas, eg    auto const V_v_def = x3::char_('v') >> '=' >> x3::int_[parser_tools::set_equals];
  auto push_back = [&](auto& ctx){ x3::_val(ctx).push_back( x3::_attr(ctx) ); };
  auto set_equals = [&](auto& ctx){ x3::_val(ctx) = x3::_attr(ctx); };

  //Use these as function objects, eg     auto const V__def = V_v[parser_tools::member_set_equals<V,std::vector<int>,&V::v>()];
  template<typename T, typename U, U T::* ptr>
  struct member_set_equals{
    template <typename Context>
    void operator()(Context const& ctx) const{
      x3::_val(ctx).*ptr = x3::_attr(ctx); 
    }
  };

  template<typename T, int reim>
  struct set_reim{
    template <typename Context>
    void operator()(Context const& ctx) const{
      reinterpret_cast<T(&)[2]>(x3::_val(ctx))[reim] = x3::_attr(ctx);
    }
  };
  
  //Error handler
  template <typename Iterator, typename Exception, typename Context>
  x3::error_handler_result
  on_error(const std::string &field, Iterator&first, Iterator const& last, Exception const& x, Context const& context){
    std::cout
      << "Error parsing " << field
      << "\nin\n" << std::string(first,last)
      << std::endl;
    return x3::error_handler_result::fail;
  }


  template<typename T>
  struct _parser_output_print: public OstreamHook{
    const T &val;
    _parser_output_print(const T&_val): val(_val){}
    void write(std::ostream &os) const{
      os << val;
    }
  };
  template<>
  void _parser_output_print<std::string>::write(std::ostream &os) const{
    os << "\"" << val << "\"";
  }
  template<>
  void _parser_output_print<bool>::write(std::ostream &os) const{
    os << (val ? "true" : "false");
  }

  
  template<typename T>
  _parser_output_print<T> parser_output_print(const T &val){
    return _parser_output_print<T>(val);
  }

  struct tabbing{
    static int & depth(){ static int d=0; return d; }
    static void increment(){ depth()++; }
    static void decrement(){ depth()--; }
    static std::string tabs(){
      std::ostringstream os;
      for(int i=0;i<depth();i++) os << '\t';
      return os.str();
    }
  };
  
};


//Bindings of types to parsers. Create overloads automatically or by hand
namespace parsers{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;

  template<typename T>
  struct parser{};

#define DEF_X3_PARSER(TYPE) \
  template<>		    \
  struct parser<TYPE>{				\
    decltype( x3::TYPE##_ ) &parse; \
  parser(): parse(x3::TYPE##_){}		\
  }
  
  DEF_X3_PARSER(int);
  DEF_X3_PARSER(double);
  DEF_X3_PARSER(bool);

#define DEF_CUSTOM_PARSER(TYPE, PARSER)		\
  template<>		    \
  struct parser<TYPE>{				\
    decltype( PARSER ) &parse; \
  parser(): parse(PARSER){}		\
  }
  
  
  //For generic vector parse we can hack the parser interface to allow templating the underlying type by static instantiating the rule definition inside the parse_rule itself.
  template<typename U>
  struct vector_T_rule{};

  template<typename T>
  struct parser< std::vector<T> >{
    x3::rule<class vector_T_rule<T>, std::vector<T> > const parse;
    parsers::parser<T> elem_parser;
    
    parser(): parse( typeid(std::vector<T>).name() ){}

    inline decltype(auto) get_def() const{
      return x3::char_('(') >> elem_parser.parse[parser_tools::push_back] >> *(',' >> elem_parser.parse[parser_tools::push_back]) > ')';
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


template<typename T>
std::ostream & operator<<(std::ostream &os, const std::vector<T> &s){
  os << '(';
  for(int i=0;i<s.size()-1;i++) os << s[i] << ", ";
  os << s.back() << ')';
  return os;
}

#include<boost/preprocessor/tuple/elem.hpp>
#include<boost/preprocessor/cat.hpp>

#define _PARSER_GEN_SEQUENCE_CREATOR_0(...)  ((__VA_ARGS__)) _PARSER_GEN_SEQUENCE_CREATOR_1
#define _PARSER_GEN_SEQUENCE_CREATOR_1(...)  ((__VA_ARGS__)) _PARSER_GEN_SEQUENCE_CREATOR_0
#define _PARSER_GEN_SEQUENCE_CREATOR_0_END
#define _PARSER_GEN_SEQUENCE_CREATOR_1_END

//Use for each for tuples
#define TUPLE_SEQUENCE_FOR_EACH(cmd, data, seq) BOOST_PP_SEQ_FOR_EACH(cmd, data , BOOST_PP_CAT(_PARSER_GEN_SEQUENCE_CREATOR_0 seq,_END) )
#define TUPLE_SEQUENCE_FOR_EACH_I(cmd, data, seq) BOOST_PP_SEQ_FOR_EACH_I(cmd, data , BOOST_PP_CAT(_PARSER_GEN_SEQUENCE_CREATOR_0 seq,_END) )
#define TUPLE_SEQUENCE_SIZE(seq) BOOST_PP_SEQ_SIZE(BOOST_PP_CAT(_PARSER_GEN_SEQUENCE_CREATOR_0 seq,_END))

#define _PARSER_MEMBER_GETNAME(elem) BOOST_PP_TUPLE_ELEM(1,elem)
#define _PARSER_MEMBER_GETTYPE(elem) BOOST_PP_TUPLE_ELEM(0,elem)
#define _PARSER_MEMBER_GETNAMESTR(elem) BOOST_PP_STRINGIZE(_PARSER_MEMBER_GETNAME(elem))

#define _PARSER_DEF_GRAMMAR_NAME(structname) BOOST_PP_CAT(structname,_grammar)

//Define parsers for the member types
#define _PARSER_MEMBER_TYPE_PARSER(elem) BOOST_PP_CAT(_PARSER_MEMBER_GETNAME(elem), _type_parse)
#define _PARSER_MEMBER_DEF_TYPE_PARSER(r,structname,elem) parsers::parser<_PARSER_MEMBER_GETTYPE(elem)> _PARSER_MEMBER_TYPE_PARSER(elem);
#define _PARSER_DEF_TYPE_PARSERS(structname, structmembers) TUPLE_SEQUENCE_FOR_EACH(_PARSER_MEMBER_DEF_TYPE_PARSER, , structmembers)

//Define and specify the rules for parsing the members
#define _PARSER_MEMBER_TAG(elem) BOOST_PP_CAT(_PARSER_MEMBER_GETNAME(elem),_parse_)
#define _PARSER_MEMBER_RULE(elem) BOOST_PP_CAT(_PARSER_MEMBER_GETNAME(elem),_parse)
#define _PARSER_MEMBER_RULESTR(elem) _PARSER_MEMBER_GETNAMESTR(elem) "_parse"
#define _PARSER_MEMBER_RULE_DEF(elem) BOOST_PP_CAT(_PARSER_MEMBER_GETNAME(elem),_parse_def)

#define _PARSER_MEMBER_DEF_RULE(r,structname,elem) \
  struct _PARSER_MEMBER_TAG(elem){					\
  template <typename Iterator, typename Exception, typename Context>	\
    x3::error_handler_result on_error(Iterator&first, Iterator const& last, Exception const& x, Context const& context){ \
    return parser_tools::on_error(BOOST_PP_STRINGIZE(structname) "." _PARSER_MEMBER_GETNAMESTR(elem),first,last,x,context); \
    }									\
  };									\
  \
  x3::rule<_PARSER_MEMBER_TAG(elem), _PARSER_MEMBER_GETTYPE(elem) > _PARSER_MEMBER_RULE(elem) = _PARSER_MEMBER_RULESTR(elem); \
  \
  auto const _PARSER_MEMBER_RULE_DEF(elem) = x3::lit(_PARSER_MEMBER_GETNAMESTR(elem)) >> '=' >> _PARSER_MEMBER_TYPE_PARSER(elem).parse[parser_tools::set_equals]; \
  \
  template <typename Iterator, typename Context, typename Attribute> \
  inline bool parse_rule( decltype(_PARSER_MEMBER_RULE(elem)) rule_ , Iterator& first, Iterator const& last , Context const& context, Attribute& attr){ \
    using boost::spirit::x3::unused; \
    static auto const def_ = (_PARSER_MEMBER_RULE(elem) = _PARSER_MEMBER_RULE_DEF(elem)); \
    return def_.parse(first, last, context, unused, attr);		\
  };

#define _PARSER_DEF_MEMBER_RULES(structname, structmembers) TUPLE_SEQUENCE_FOR_EACH(_PARSER_MEMBER_DEF_RULE, structname, structmembers)

//Define the rule for the main structure
#define _PARSER_STRUCT_RULE_TAG(name) BOOST_PP_CAT(name,_)
#define _PARSER_STRUCT_RULE(name) BOOST_PP_CAT(name,_)
#define _PARSER_STRUCT_RULE_DEF(name) BOOST_PP_CAT(name,__def)

#define _PARSER_DEF_STRUCT_RULE_MEMBER_GEN(r,structname,elem) \
  >> _PARSER_MEMBER_RULE(elem)[parser_tools::member_set_equals<structname,_PARSER_MEMBER_GETTYPE(elem),& structname :: _PARSER_MEMBER_GETNAME(elem)>()]

#define _PARSER_DEF_STRUCT_RULE_DEF(NAME)\
  struct _PARSER_STRUCT_RULE_TAG(NAME){\
    template <typename Iterator, typename Exception, typename Context>	\
      x3::error_handler_result on_error(Iterator&first, Iterator const& last, Exception const& x, Context const& context){ \
      return parser_tools::on_error(BOOST_PP_STRINGIZE(NAME),first,last,x,context); \
    }									\
  };									\
									\
									\
  x3::rule<_PARSER_STRUCT_RULE_TAG(NAME), NAME> const _PARSER_STRUCT_RULE(NAME) = BOOST_PP_STRINGIZE(_PARSER_STRUCT_RULE(NAME));

//Specify the rules for parsing the structure
#define _PARSER_DEF_STRUCT_RULE_IMPL(NAME,MEMSEQ) \
  auto const _PARSER_STRUCT_RULE_DEF(NAME) = x3::eps > '{'		\
  TUPLE_SEQUENCE_FOR_EACH(_PARSER_DEF_STRUCT_RULE_MEMBER_GEN, NAME, MEMSEQ) \
    > '}'; \
  BOOST_SPIRIT_DEFINE(_PARSER_STRUCT_RULE(NAME));

//Register the parser
#define _PARSER_DEF_ADD_PARSER_TO_NAMESPACE(NAME) \
  namespace parsers{ \
  template<> \
  struct parser<NAME>{\
    decltype( _PARSER_DEF_GRAMMAR_NAME(NAME)::_PARSER_STRUCT_RULE(NAME) ) &parse; \
    parser(): parse(_PARSER_DEF_GRAMMAR_NAME(NAME)::_PARSER_STRUCT_RULE(NAME) ){} \
  }; \
};

//Write operator<< and operator>>
#define _PARSER_DEF_OSTREAM_MEMBER_WRITE(r,structname,elem)  os << parser_tools::tabbing::tabs() << BOOST_PP_STRINGIZE(_PARSER_MEMBER_GETNAME(elem)) " = " << parser_tools::parser_output_print(s._PARSER_MEMBER_GETNAME(elem)) << std::endl;

#define _PARSER_DEF_DEFINE_OSTREAM_WRITE(NAME, MEMSEQ) \
  std::ostream & operator<<(std::ostream &os, const NAME &s){ \
  os << "{\n"; parser_tools::tabbing::increment(); \
  TUPLE_SEQUENCE_FOR_EACH(_PARSER_DEF_OSTREAM_MEMBER_WRITE, NAME, MEMSEQ) \
  parser_tools::tabbing::decrement(); os << parser_tools::tabbing::tabs() << "}";  \
  return os; \
}

#define _PARSER_DEF_DEFINE_ISTREAM_READ(NAME, MEMSEQ) \
std::istream & operator>>(std::istream &is, NAME &s){ \
  namespace ascii = boost::spirit::x3::ascii; \
  namespace x3 = boost::spirit::x3; \
  \
  is >> std::noskipws; \
  boost::spirit::istream_iterator f(is), l; \
  parsers::parser<NAME> vp; \
  \
  bool r = x3::phrase_parse(f, l, vp.parse, ascii::space, s); \
  \
  if(!r){ \
    throw std::runtime_error("Parsing of type " BOOST_PP_STRINGIZE(NAME) " failed\n"); \
  } \
  return is; \
}

//Combine all the above the specify a parser
//To use, define a structure as
// GENERATE_PARSER( <STRUCT_NAME> ,
// 		 (TYPE1, MEMBER1)
// 		 (TYPE2, MEMBER2).... )

#define GENERATE_PARSER(structname, structmembers)					 \
  namespace _PARSER_DEF_GRAMMAR_NAME(structname){				\
    namespace ascii = boost::spirit::x3::ascii; \
    namespace x3 = boost::spirit::x3; \
    _PARSER_DEF_TYPE_PARSERS(structname, structmembers) \
    _PARSER_DEF_MEMBER_RULES(structname, structmembers) \
   \
   _PARSER_DEF_STRUCT_RULE_DEF(structname) \
   \
   _PARSER_DEF_STRUCT_RULE_IMPL(structname,structmembers) \
  }; \
  _PARSER_DEF_ADD_PARSER_TO_NAMESPACE(structname) \
  _PARSER_DEF_DEFINE_OSTREAM_WRITE(structname,structmembers) \
  _PARSER_DEF_DEFINE_ISTREAM_READ(structname,structmembers)


//Convenience function for generating members using the same struct as the parser
#define _PARSER_DEF_MEMBER_DEFINE(r,data,elem) _PARSER_MEMBER_GETTYPE(elem) _PARSER_MEMBER_GETNAME(elem);
#define GENERATE_MEMBERS(structmembers) TUPLE_SEQUENCE_FOR_EACH(_PARSER_DEF_MEMBER_DEFINE,   ,structmembers)




//Some useful macros

//For a sequence of tuples, output a new sequence of tuples with element 'idx' removed 
// #define _POPI(r,idx,elem) BOOST_PP_TUPLE_REMOVE(elem,idx)
// #define SEQPOP(seq,idx) TUPLE_SEQUENCE_FOR_EACH(

// BOOST_PP_SEQ_FOR_EACH(_POP1, idx  , BOOST_PP_CAT(_PARSER_GEN_SEQUENCE_CREATOR_0 seq,_END) )





//Example hand-written parser for simple struct
/*
struct V{
  std::vector<int> v;
};
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


#endif
