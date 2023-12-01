#ifndef _JACKKNIFE_DIST_ITERATE_H_
#define _JACKKNIFE_DIST_ITERATE_H_

#include<config.h>
#include<distribution/jackknife/class.h>
#include<distribution/distribution_iterate.h>

SARLAC_START_NAMESPACE

template<typename T, template<typename> class V>
struct iterate<jackknifeDistribution<T,V> >{
  typedef T type;

  static inline int size(const jackknifeDistribution<T,V> &from){ return from.size(); } 
  static inline const T& at(const int i, const jackknifeDistribution<T,V> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, jackknifeDistribution<T,V> &from){
    return from.sample(i);
  }
  static inline std::vector<int> unmap(const int i, const jackknifeDistribution<T,V> &from){ 
    return std::vector<int>({i});
  }
};

SARLAC_END_NAMESPACE
#endif
