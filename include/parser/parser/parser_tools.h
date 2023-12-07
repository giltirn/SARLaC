#ifndef _SARLAC_PARSER_TOOLS_H_
#define _SARLAC_PARSER_TOOLS_H_

#include <vector>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<parser/parser/parser_macros.h>

SARLAC_START_NAMESPACE

namespace parser_tools{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  
  //Use these as lambdas, eg    auto const V_v_def = x3::char_('v') >> '=' >> x3::int_[parser_tools::set_equals];
  inline auto push_back = [](auto& ctx){ x3::_val(ctx).push_back( x3::_attr(ctx) ); };
  inline auto set_equals = [](auto& ctx){ x3::_val(ctx) = x3::_attr(ctx); };

  //Use these as function objects, eg     auto const V__def = V_v[parser_tools::member_set_equals<V,std::vector<int>,&V::v>()];
  template<typename T, typename U, U T::* ptr>
  struct member_set_equals{
    template <typename Context>
    inline void operator()(Context const& ctx) const{
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
  
  template<typename T>
  struct set_elem{
    int elem;
    set_elem(int elem): elem(elem){}

    template <typename Context>
    inline void operator()(Context const& ctx) const{
      x3::_val(ctx)[elem] = x3::_attr(ctx);
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

  //Tools for writing various types in the format used by the parser
  template<typename T>
  struct _parser_output_print: public OstreamHook{
    const T &val;
    _parser_output_print(const T&_val): val(_val){}
    void write(std::ostream &os) const{
      os << val;
    }
  };

  template<typename T>
  inline _parser_output_print<T> parser_output_print(const T &val){
    return _parser_output_print<T>(val);
  }
  template<typename T>
  inline void parser_output_print(std::ostream &os, const T &val){
    _parser_output_print<T> p(val); p.write(os);
  }

  template<>
  inline void _parser_output_print<std::string>::write(std::ostream &os) const{
    os << "\"" << val << "\"";
  }
  template<>
  inline void _parser_output_print<bool>::write(std::ostream &os) const{
    os << (val ? "true" : "false");
  }

  //std::vector
  template<typename T>
  struct _parser_output_print< std::vector<T> >: public OstreamHook{
    const std::vector<T> &val;
    _parser_output_print(const std::vector<T> &_val): val(_val){}
    void write(std::ostream &os) const{
      os << '(';
      if(val.size() > 0){
	for(int i=0;i<val.size()-1;i++){
	  parser_output_print(os, val[i]);
	  os << ", ";
	}
	parser_output_print(os, val.back());
      }
      os << ')';
    }
  };

  //std::array
  template<typename T, size_t N>
  struct _parser_output_print< std::array<T,N> >: public OstreamHook{
    const std::array<T,N> &val;
    _parser_output_print(const std::array<T,N> &_val): val(_val){}
    void write(std::ostream &os) const{
      os << '(';
      if(N > 0){
	for(int i=0;i<N-1;i++){
	  parser_output_print(os, val[i]);
	  os << ", ";
	}
	parser_output_print(os, val.back());
      }
      os << ')';
    }
  };

  //std::pair
  template<typename T, typename U>
  struct _parser_output_print< std::pair<T,U> >: public OstreamHook{
    const std::pair<T,U> &val;
    _parser_output_print(const std::pair<T,U> &_val): val(_val){}
    void write(std::ostream &os) const{
      os << '{';
      parser_output_print(os, val.first);
      os << ", ";
      parser_output_print(os, val.second);
      os << "}";
    }
  };

  //std::map
  template<typename T, typename N>
  struct _parser_output_print< std::map<T,N> >: public OstreamHook{
    const std::map<T,N> &val;
    bool in_line;
    _parser_output_print(const std::map<T,N> &_val, const bool in_line = false): val(_val), in_line(in_line){}
    void write(std::ostream &os) const{
      auto last_it = std::next(val.begin(), val.size()-1);
      os << "{"; if(!in_line) os << std::endl;
      for(auto it = val.begin(); it != val.end(); it++){
	parser_output_print(os, it->first);
	os << " : ";
	parser_output_print(os, it->second);
	if(it != last_it) os << ",";
	if(!in_line) os << std::endl;
      }
      os << "}";
    }
  };


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

SARLAC_END_NAMESPACE
#endif
