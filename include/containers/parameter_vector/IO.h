#ifndef _CPSFIT_PARAMETER_VECTOR_IO_H_
#define _CPSFIT_PARAMETER_VECTOR_IO_H_

#include<config.h>
#include<utils/macros.h>
#include<parser/parser.h>
#include<containers/parameter_vector/class.h>

CPSFIT_START_NAMESPACE

#ifdef HAVE_HDF5
template<typename D>
inline void write(HDF5writer &writer, const parameterVector<D> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, parameterVector<D> &d, const std::string &tag){ d.read(reader,tag); }
#endif


/* template<typename Numeric>  */
/* std::ostream & operator<<(std::ostream & stream, const parameterVector<Numeric> &vec){ */
/*   stream << "("; */
/*   for(int i=0;i<vec.size();i++) */
/*     stream << vec[i] << (i != vec.size()-1 ? " " : ")"); */
/*   return stream; */
/* } */


//Define the X3 parser
namespace parsers{
  template<typename U>
  struct parameter_vector_T_rule{};

  template<typename T>
  struct parser< parameterVector<T> >{
    x3::rule<class parameter_vector_T_rule<T>, parameterVector<T> > const parse;
    parsers::parser<T> elem_parser;
    
    parser(): parse( typeid(parameterVector<T>).name() ){}

    inline decltype(auto) get_def() const{
      return x3::char_('(') >> *( elem_parser.parse[parser_tools::push_back] >> *(',' >> elem_parser.parse[parser_tools::push_back]) ) > ')';
    }
  };
  template <typename T, typename Iterator, typename Context, typename Attribute>
  bool parse_rule(x3::rule<class parameter_vector_T_rule<T>, parameterVector<T> > const rule_ , Iterator& first, Iterator const& last , Context const& context, Attribute& attr){
    using boost::spirit::x3::unused;
    static parser<parameterVector<T> > inst;
    static auto const parse_def = inst.get_def();
    static auto const def_ = (inst.parse = parse_def);
    return def_.parse(first, last, context, unused, attr);
  }
};

CPSFIT_END_NAMESPACE

#endif
