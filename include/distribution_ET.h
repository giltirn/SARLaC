#ifndef _DISTRIBUTION_ET_H
#define _DISTRIBUTION_ET_H

#include<cstdio>
#include<vector>
#include<type_traits>
#include<cassert>
#include<generic_ET.h>

template<typename A>
struct getElem<distribution<A> >{
  static inline auto elem(const distribution<A> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(distribution<A> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const distribution<A> &v){ return v.size(); }
};

template<typename A>
struct getElem<jackknifeDistribution<A> >{
  static inline auto elem(const jackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(jackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const jackknifeDistribution<A> &v){ return v.size(); }
};

template<typename A>
struct getElem<doubleJackknifeDistribution<A> >{
  static inline auto elem(const doubleJackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(doubleJackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const doubleJackknifeDistribution<A> &v){ return v.size(); }
};


#endif
