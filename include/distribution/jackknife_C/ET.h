#ifndef _JACKKNIFE_C_DIST_ET_H_
#define _JACKKNIFE_C_DIST_ET_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/jackknife_C/class.h>

SARLAC_START_NAMESPACE

template<typename A, template<typename> class V>
struct getElem<jackknifeCdistribution<A,V> >{
  static inline auto elem(const jackknifeCdistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(jackknifeCdistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const jackknifeCdistribution<A,V> &v){ return v.size(); }
};

SARLAC_END_NAMESPACE
#endif
