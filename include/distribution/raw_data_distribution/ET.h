#ifndef _RAW_DATA_DIST_ET_H_
#define _RAW_DATA_DIST_ET_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/raw_data_distribution/class.h>

SARLAC_START_NAMESPACE

template<typename A, template<typename> class V>
  struct getElem<rawDataDistribution<A,V> >{
  static inline auto elem(const rawDataDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(rawDataDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline size_t common_properties(const rawDataDistribution<A,V> &v){ return v.size(); }
};

SARLAC_END_NAMESPACE
#endif
