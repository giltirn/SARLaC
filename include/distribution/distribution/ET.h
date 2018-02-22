#ifndef _DISTRIBUTION_ET_H_
#define _DISTRIBUTION_ET_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/distribution/class.h>

CPSFIT_START_NAMESPACE

template<typename A, template<typename> class V>
  struct getElem<distribution<A,V> >{
  static inline auto elem(const distribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(distribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const distribution<A,V> &v){ return v.size(); }
};

CPSFIT_END_NAMESPACE
#endif
