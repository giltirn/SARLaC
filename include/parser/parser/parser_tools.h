#ifndef _CPSFIT_PARSER_TOOLS_H_
#define _CPSFIT_PARSER_TOOLS_H_

#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<parser/parser/parser_macros.h>

CPSFIT_START_NAMESPACE

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
    // std::cout
    //   << "Error parsing " << field
    //   << "\nin\n" << std::string(first,last)
    //   << std::endl;
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

CPSFIT_END_NAMESPACE
#endif
