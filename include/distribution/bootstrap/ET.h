#ifndef _BOOTSTRAP_DIST_ET_H_
#define _BOOTSTRAP_DIST_ET_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/bootstrap/class.h>

SARLAC_START_NAMESPACE

template<typename A, template<typename> class V>
struct getElem<bootstrapDistribution<A,V> >{
  static inline auto elem(const bootstrapDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return i==-1 ? v.best() : v.sample(i); }
  static inline auto elem(bootstrapDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return i==-1 ? v.best() : v.sample(i); }
  static inline typename bootstrapDistribution<A,V>::initType common_properties(const bootstrapDistribution<A,V> &v){ return typename bootstrapDistribution<A,V>::initType(v.size(),v.confidence()); }
};

SARLAC_END_NAMESPACE
#endif
