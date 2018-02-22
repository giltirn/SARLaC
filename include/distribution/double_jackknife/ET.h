#ifndef _DOUBLE_JACKKNIFE_DIST_ET_H_
#define _DOUBLE_JACKKNIFE_DIST_ET_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/double_jackknife/class.h>

CPSFIT_START_NAMESPACE

template<typename A, template<typename> class V>
  struct getElem<doubleJackknifeDistribution<A,V> >{
  static inline auto elem(const doubleJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(doubleJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const doubleJackknifeDistribution<A,V> &v){ return v.size(); }
};

CPSFIT_END_NAMESPACE
#endif
