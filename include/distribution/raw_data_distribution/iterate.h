#ifndef _RAW_DATA_DIST_ITERATE_H_
#define _RAW_DATA_DIST_ITERATE_H_

#include<config.h>
#include<distribution/raw_data_distribution/class.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE

template<typename T, template<typename> class V>
struct iterate<rawDataDistribution<T,V> >{
  typedef T type;

  static inline int size(const rawDataDistribution<T,V> &from){ return from.size(); } 
  static inline const T& at(const int i, const rawDataDistribution<T,V> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, rawDataDistribution<T,V> &from){
    return from.sample(i);
  }
  static inline std::vector<int> unmap(const int i, const rawDataDistribution<T,V> &from){ 
    return std::vector<int>({i});
  }
};

CPSFIT_END_NAMESPACE
#endif
