#ifndef _SARLAC_PARSER_MACROS_H
#define _SARLAC_PARSER_MACROS_H

//Here be dragons!
/*                                 _/|__ */
/*             _,-------,        _/ -|  \_     /~>. */
/*          _-~ __--~~/\ |      (  \   /  )   | / | */
/*       _-~__--    //   \\      \ *   * /   / | || */
/*    _-~_--       //     ||      \     /   | /  /| */
/*   ~ ~~~~-_     //       \\     |( " )|  / | || / */
/*           \   //         ||    | VWV | | /  /// */
/*     |\     | //           \\ _/      |/ | ./| */
/*     | |    |// __         _-~         \// |  / */
/*    /  /   //_-~  ~~--_ _-~  /          |\// / */
/*   |  |   /-~        _-~    (     /   |/ / / */
/*  /   /           _-~  __    |   |____|/ */
/* |   |__         / _-~  ~-_  (_______  `\ */
/* |      ~~--__--~ /  _     \        __\))) */
/*  \               _-~       |     ./  \ */
/*   ~~--__        /         /    _/     | */
/*         ~~--___/       _-_____/      / */
/*          _____/     _-_____/      _-~ */
/*       /^<  ___       -____         -____ */
/*          ~~   ~~--__      ``\--__       ``\ */
/*                     ~~--\)\)\)   ~~--\)\)\) */
//....
//Specifically utility macros for generating enum and struct parsers using boost spirit X3
//Users should not invoke these directly

#include<boost/preprocessor/tuple/elem.hpp>
#include<boost/preprocessor/cat.hpp>

#include<config.h>
#include<utils/macros.h>
#include<parser/parser/parsers.h>

#define DEF_X3_PARSER(TYPE) \
  template<>		    \
  struct parser<TYPE>{				\
    decltype( x3::TYPE##_ ) &parse; \
  parser(): parse(x3::TYPE##_){}		\
  }

#define DEF_CUSTOM_PARSER(TYPE, PARSER)		\
  template<>		    \
  struct parser<TYPE>{				\
    decltype( PARSER ) &parse; \
  parser(): parse(PARSER){}		\
  }

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
#define _PARSER_MEMBER_DEF_TYPE_PARSER(r,dummy,elem) ::SARLaC::parsers::parser<_PARSER_MEMBER_GETTYPE(elem)> _PARSER_MEMBER_TYPE_PARSER(elem);
#define _PARSER_DEF_TYPE_PARSERS(structmembers) TUPLE_SEQUENCE_FOR_EACH(_PARSER_MEMBER_DEF_TYPE_PARSER, , structmembers)

#define _PARSER_MEMBER_TYPE_PARSER_INST(elem) ::SARLaC::parsers::parser_instance<_PARSER_MEMBER_GETTYPE(elem)>::get() 


//Define and specify the rules for parsing the members
#define _PARSER_MEMBER_TAG(elem) BOOST_PP_CAT(_PARSER_MEMBER_GETNAME(elem),_parse_)
#define _PARSER_MEMBER_RULE(elem) BOOST_PP_CAT(_PARSER_MEMBER_GETNAME(elem),_parse)
#define _PARSER_MEMBER_RULESTR(elem) _PARSER_MEMBER_GETNAMESTR(elem) "_parse"
#define _PARSER_MEMBER_RULE_DEF(elem) BOOST_PP_CAT(_PARSER_MEMBER_GETNAME(elem),_parse_def)

#define _PARSER_MEMBER_DEF_RULE(r,structname,elem) \
  struct _PARSER_MEMBER_TAG(elem): ::SARLaC::parsers::error_handler{ }; \
  \
  inline auto & _PARSER_MEMBER_RULE(elem)(){				\
     static x3::rule<_PARSER_MEMBER_TAG(elem), _PARSER_MEMBER_GETTYPE(elem) > rule_inst = _PARSER_MEMBER_RULESTR(elem); \
     return rule_inst; \
  }\
  template <typename Iterator, typename Context, typename Attribute> \
  inline bool parse_rule( x3::rule<_PARSER_MEMBER_TAG(elem), _PARSER_MEMBER_GETTYPE(elem) > rule_ , Iterator& first, Iterator const& last , Context const& context, Attribute& attr){ \
    using boost::spirit::x3::unused; \
    static auto const _PARSER_MEMBER_RULE_DEF(elem) = x3::lit(_PARSER_MEMBER_GETNAMESTR(elem)) > '=' > _PARSER_MEMBER_TYPE_PARSER_INST(elem).parse[::SARLaC::parser_tools::set_equals]; \
    static auto const def_ = (_PARSER_MEMBER_RULE(elem)() = _PARSER_MEMBER_RULE_DEF(elem)); \
    return def_.parse(first, last, context, unused, attr);		\
  };

#define _PARSER_DEF_MEMBER_RULES(structname, structmembers) TUPLE_SEQUENCE_FOR_EACH(_PARSER_MEMBER_DEF_RULE, structname, structmembers)

//Define the rule for the main structure
#define _PARSER_DEF_STRUCT_RULE_MEMBER_GEN(r,structname,elem) \
  > _PARSER_MEMBER_RULE(elem)()[::SARLaC::parser_tools::member_set_equals<structname,_PARSER_MEMBER_GETTYPE(elem),& structname :: _PARSER_MEMBER_GETNAME(elem)>()]

#define _PARSER_DEF_STRUCT_RULE_DEF(NAME)\
  struct main_rule_handler: ::SARLaC::parsers::error_handler{ };	\
  									\
  inline auto & main_rule(){						\
     static x3::rule<main_rule_handler, NAME> const main_rule = BOOST_PP_STRINGIZE(NAME); \
     return main_rule; \
  }

//Specify the rules for parsing the structure
#define _PARSER_DEF_STRUCT_RULE_IMPL(NAME,MEMSEQ) \
  template <typename Iterator, typename Context, typename Attribute> \
  inline bool parse_rule( x3::rule<main_rule_handler, NAME> rule_ , Iterator& first, Iterator const& last , Context const& context, Attribute& attr){ \
    using boost::spirit::x3::unused;					\
    static auto const main_rule_def = x3::char_('{')			\
      TUPLE_SEQUENCE_FOR_EACH(_PARSER_DEF_STRUCT_RULE_MEMBER_GEN, NAME, MEMSEQ) \
      > '}';								\
    static auto const def_ = (main_rule() = main_rule_def);		\
    return def_.parse(first, last, context, unused, attr);		\
  };


//Register the parser
#define _PARSER_DEF_ADD_PARSER_TO_NAMESPACE(NAME,GRAMMAR)	\
    template<>							\
    struct SARLaC::parsers::parser<NAME>{				\
      const typename std::decay<decltype(GRAMMAR::main_rule())>::type &parse;	\
      parser(): parse(GRAMMAR::main_rule() ){}			\
    };

//Write operator<< and operator>>
#define _PARSER_DEF_OSTREAM_MEMBER_WRITE(r,structname,elem)  os << ::SARLaC::parser_tools::tabbing::tabs() << BOOST_PP_STRINGIZE(_PARSER_MEMBER_GETNAME(elem)) " = " << ::SARLaC::parser_tools::parser_output_print(s._PARSER_MEMBER_GETNAME(elem)) << std::endl;

#define _PARSER_DEF_DEFINE_OSTREAM_WRITE(NAME, MEMSEQ) \
  inline std::ostream & operator<<(std::ostream &os, const NAME &s){	\
    os << "{\n"; ::SARLaC::parser_tools::tabbing::increment();		\
    TUPLE_SEQUENCE_FOR_EACH(_PARSER_DEF_OSTREAM_MEMBER_WRITE, NAME, MEMSEQ) \
      ::SARLaC::parser_tools::tabbing::decrement(); os << ::SARLaC::parser_tools::tabbing::tabs() << "}"; \
    return os;								\
  }

#define _GENERATE_PARSER_GRAMMAR(structname, grammar, structmembers)		\
    namespace grammar{							\
      namespace ascii = boost::spirit::x3::ascii;			\
      namespace x3 = boost::spirit::x3;					\
      _PARSER_DEF_MEMBER_RULES(structname, structmembers)		\
									\
      _PARSER_DEF_STRUCT_RULE_DEF(structname)				\
									\
      _PARSER_DEF_STRUCT_RULE_IMPL(structname,structmembers)		\
    };	

//Parser generated from outside SARLaC namespace
#define _GENERATE_PARSER_GM(structname, grammar, structmembers)		\
    _GENERATE_PARSER_GRAMMAR(structname, grammar, structmembers)	\
    _PARSER_DEF_ADD_PARSER_TO_NAMESPACE(structname,grammar)		\
    _PARSER_DEF_DEFINE_OSTREAM_WRITE(structname,structmembers)

//Generate the parser from outside SARLaC namespace
#define _GENERATE_PARSER(structname, structmembers) _GENERATE_PARSER_GM(structname, _PARSER_DEF_GRAMMAR_NAME(structname), structmembers)


//Convenience function for generating members using the same struct as the parser
#define _PARSER_DEF_MEMBER_DEFINE(r,data,elem) _PARSER_MEMBER_GETTYPE(elem) _PARSER_MEMBER_GETNAME(elem);
#define _GENERATE_MEMBERS(structmembers) TUPLE_SEQUENCE_FOR_EACH(_PARSER_DEF_MEMBER_DEFINE,   ,structmembers)


//Similar macros for generating parser for enums
#define _GEN_ENUM_ENUMDEF(enumname, enummembers) enum class enumname { BOOST_PP_SEQ_ENUM(enummembers) };

#define _GEN_ENUM_STR_E(r,data,i,elem) BOOST_PP_COMMA_IF(i) BOOST_PP_STRINGIZE(elem)
#define _GEN_ENUM_STR(enumname, enummembers) \
  inline std::string toString(const enumname d){				\
    const static std::vector<std::string> str = { BOOST_PP_SEQ_FOR_EACH_I(_GEN_ENUM_STR_E,  ,  enummembers) }; \
    int dd = static_cast<int>(d);					\
    if(dd < 0 || dd >= str.size()){ \
      std::ostringstream os; os << "Unknown " BOOST_PP_STRINGIZE(enumname) " idx " << dd; \
      return os.str(); \
    }else return str[static_cast<int>(d)]; \
  }

#define _GEN_ENUM_PARSER_MATCH(enumname, enummembers) \
  struct BOOST_PP_CAT(enumname,_match){ \
    template <typename Context> \
    void operator()(Context const& ctx) const{ \
      const static std::vector<std::string> str = { BOOST_PP_SEQ_FOR_EACH_I(_GEN_ENUM_STR_E,  ,  enummembers) }; \
      std::string tag = x3::_attr(ctx); \
      enumname &val = x3::_val(ctx);	\
      for(int i=0;i<str.size();i++) \
	if(tag == str[i]){ val = static_cast<enumname>(i); return; }	\
      ::SARLaC::error_exit(std::cout << "Unknown " BOOST_PP_STRINGIZE(enumname) " : " << tag << std::endl); \
    } \
  };

#define _GEN_ENUM_PARSER_BODY(enumname, enummembers) \
  namespace BOOST_PP_CAT(enumname,_parser){	     \
    namespace ascii = boost::spirit::x3::ascii; \
    namespace x3 = boost::spirit::x3;\
    auto const enumparse = x3::lexeme[+x3::char_("a-zA-Z0-9_")]; \
    struct BOOST_PP_CAT(enumname,_tag) : ::SARLaC::parsers::error_handler{ }; \
    _GEN_ENUM_PARSER_MATCH(enumname, enummembers) \
    x3::rule<BOOST_PP_CAT(enumname,_tag), enumname> const BOOST_PP_CAT(enumname,_) = BOOST_PP_STRINGIZE(BOOST_PP_CAT(enumname,_)); \
    auto const BOOST_PP_CAT(enumname, __def) = enumparse[BOOST_PP_CAT(enumname,_match)()]; \
    BOOST_SPIRIT_DEFINE(BOOST_PP_CAT(enumname,_)); \
  }; \
    template<>		   \
    struct SARLaC::parsers::parser<enumname>{				\
      decltype( BOOST_PP_CAT(enumname, _parser)::BOOST_PP_CAT(enumname, _) ) &parse; \
      parser(): parse( BOOST_PP_CAT(enumname, _parser)::BOOST_PP_CAT(enumname, _) ){} \
    }; \
    

#define _GEN_ENUM_DEFINE_OSTREAM_WRITE(enumname, enummembers)	      \
  std::ostream & operator<<(std::ostream &os, const enumname s){ \
    os << toString(s); \
    return os; \
  }



//Generate the parser from outside SARLaC namespace
#define _GENERATE_ENUM_PARSER(enumname, enummembers)	\
  _GEN_ENUM_STR(enumname, enummembers)				\
  _GEN_ENUM_PARSER_BODY(enumname, enummembers)			\
  _GEN_ENUM_DEFINE_OSTREAM_WRITE(enumname, enummembers)		\

//Define an enum and a parser to go along with it (from outside SARLaC namespace)
//To use  GENERATE_ENUM_AND_PARSER( <ENUM NAME>, (<ELEM1>)(<ELEM2>)(<ELEM3>)... )
#define _GENERATE_ENUM_AND_PARSER(enumname, enummembers)	\
  _GEN_ENUM_ENUMDEF(enumname, enummembers) \
  _GENERATE_ENUM_PARSER(enumname, enummembers)
  
#endif
