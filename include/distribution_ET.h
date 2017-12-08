#ifndef _DISTRIBUTION_ET_H
#define _DISTRIBUTION_ET_H

#include<cstdio>
#include<vector>
#include<type_traits>
#include<cassert>
#include<generic_ET.h>
#include<distribution.h>

template<typename A, template<typename> class V>
  struct getElem<distribution<A,V> >{
  static inline auto elem(const distribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(distribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const distribution<A,V> &v){ return v.size(); }
};

template<typename A, template<typename> class V>
  struct getElem<rawDataDistribution<A,V> >{
  static inline auto elem(const rawDataDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(rawDataDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const rawDataDistribution<A,V> &v){ return v.size(); }
};

template<typename A, template<typename> class V>
  struct getElem<jackknifeDistribution<A,V> >{
  static inline auto elem(const jackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(jackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const jackknifeDistribution<A,V> &v){ return v.size(); }
};

template<typename A, template<typename> class V>
  struct getElem<jackknifeCdistribution<A,V> >{
  static inline auto elem(const jackknifeCdistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(jackknifeCdistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const jackknifeCdistribution<A,V> &v){ return v.size(); }
};


template<typename A, template<typename> class V>
  struct getElem<doubleJackknifeDistribution<A,V> >{
  static inline auto elem(const doubleJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(doubleJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const doubleJackknifeDistribution<A,V> &v){ return v.size(); }
};


#endif
