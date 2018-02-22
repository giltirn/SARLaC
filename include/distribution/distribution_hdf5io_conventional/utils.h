#ifndef _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_UTILS_H_
#define _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_UTILS_H_

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

template<typename DistributionType>
inline std::string getDistributionTypeString(){
  std::string base = printType<DistributionType>();
  size_t p = base.find_first_of(','); //ignore vector type
  if(p == std::string::npos) return base;
  else{
    std::ostringstream os; os << base.substr(0, p) << '>';
    return os.str();
  }
}

CPSFIT_END_NAMESPACE
#endif
