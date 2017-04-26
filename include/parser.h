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
  
  //Rule for std::string
  class string_rule_{
    template <typename Iterator, typename Exception, typename Context>
    x3::error_handler_result on_error(Iterator&first, Iterator const& last, Exception const& x, Context const& context){
      return parser_tools::on_error("string_rule",first,last,x,context);
    }
  };
  
  auto const quoted_string = x3::lexeme['"' > *(x3::char_ - '"') > '"'];
  x3::rule<string_rule_, std::string> const string_rule = "string_rule";
  auto const string_rule_def = quoted_string[parser_tools::set_equals];
  BOOST_SPIRIT_DEFINE(string_rule);
  
  template<>
  struct parser<std::string>{
    decltype( string_rule ) &parse;
    parser(): parse(string_rule){}
  };

};


template<typename T>
std::ostream & operator<<(std::ostream &os, const std::vector<T> &s){
  os << '(';
  for(int i=0;i<s.size()-1;i++) os << s[i] << ", ";
  os << s.back() << ')';
  return os;
}


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
