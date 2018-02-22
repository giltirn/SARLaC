#ifndef _XML_SERIALIZE_XML_TRANSLATOR_H_
#define _XML_SERIALIZE_XML_TRANSLATOR_H_

#include<vector>
#include<sstream>
#include<boost/optional.hpp>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//Space-separated POD array
template<typename PODtype, typename std::enable_if<std::is_pod<PODtype>::value, int>::type = 0>
struct pod_vector_translator
{
  typedef std::string internal_type;
  typedef std::vector<PODtype> external_type;
  boost::optional<external_type> get_value(const internal_type& str){
    boost::optional<external_type> out(true,external_type());
    std::stringstream ss(str);
    PODtype tmp;
    while(ss){
      if(ss.eof()) break;
      ss >> tmp;
      out.get().push_back(tmp);
    }
    return out;
  }
};


CPSFIT_END_NAMESPACE


#endif
