#ifndef _SUPERMULTI_DIST_ITERATE_H_
#define _SUPERMULTI_DIST_ITERATE_H_

#include<config.h>
#include<distribution/supermulti/class.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE

template<typename T>
struct iterate<superMultiDistribution<T> >{
  typedef T type;

  static inline int size(const superMultiDistribution<T> &from){ return from.size()+1; } 
  static inline const T& at(const int i, const superMultiDistribution<T> &from){
    return i==0 ? from.best() : from.sample(i-1);
  }
  static inline T & at(const int i, superMultiDistribution<T> &from){
    return i==0 ? from.best() : from.sample(i-1);
  }
  static inline std::vector<int> unmap(const int i, const superMultiDistribution<T> &from){ 
    return std::vector<int>({i});
  }
  
};

CPSFIT_END_NAMESPACE
#endif
