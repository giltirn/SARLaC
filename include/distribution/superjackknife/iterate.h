#ifndef _SUPER_JACKKNIFE_DIST_ITERATE_H_
#define _SUPER_JACKKNIFE_DIST_ITERATE_H_

#include<config.h>
#include<distribution/superjackknife/class.h>
#include<distribution/distribution_iterate.h>

SARLAC_START_NAMESPACE

template<typename T>
struct iterate<superJackknifeDistribution<T> >{
  typedef T type;

  static inline int size(const superJackknifeDistribution<T> &from){ return from.size()+1; } 
  static inline const T& at(const int i, const superJackknifeDistribution<T> &from){
    return i==0 ? from.best() : from.sample(i-1);
  }
  static inline T & at(const int i, superJackknifeDistribution<T> &from){
    return i==0 ? from.best() : from.sample(i-1);
  }
  static inline std::vector<int> unmap(const int i, const superJackknifeDistribution<T> &from){ 
    return std::vector<int>({i});
  }
  
};

SARLAC_END_NAMESPACE
#endif
