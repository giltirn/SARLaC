#ifndef _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_UTILS_H_
#define _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_UTILS_H_

#include<regex>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<utils/template_wizardry.h>

CPSFIT_START_NAMESPACE

template<typename T, int hasResize>
struct _actionResize{
  static inline void doit(T &v, const int n){}
};

template<typename T>
struct _actionResize<T,1>{
  static inline void doit(T &v, const int n){
    v.resize(n);
  }
};

//Stringify the distrubution such that   distributionType<T,vector_type>  becomes  "distributionType<T>"
//The vector type is irrelevant for the purposes of saving and loading distributions
template<typename DistributionType>
inline std::string getDistributionTypeString(){
  std::string tpstr = printType<DistributionType>();
  
  std::regex r(R"(^(?:\w+\:\:)?(\w+)\<(\w+),)", std::regex_constants::ECMAScript);
  std::smatch m;
  if(std::regex_search(tpstr, m, r)){
    std::ostringstream os;
    os << m[1].str() << "<" << m[2].str() << ">";
    return os.str();
  }else{
    error_exit(std::cout << "getDistributionTypeString could not parse type " << tpstr << std::endl);
  }
}

CPSFIT_END_NAMESPACE
#endif
