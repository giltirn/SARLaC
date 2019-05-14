#ifndef _BLOCK_DOUBLE_JACKKNIFE_DIST_ET_H_
#define _BLOCK_DOUBLE_JACKKNIFE_DIST_ET_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/block_double_jackknife/class.h>

CPSFIT_START_NAMESPACE

template<typename A, template<typename> class V>
  struct getElem<blockDoubleJackknifeDistribution<A,V> >{
  static inline auto elem(const blockDoubleJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(blockDoubleJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline std::pair<int,int> common_properties(const blockDoubleJackknifeDistribution<A,V> &v){ return std::pair<int,int>(v.nSamplesUnbinned(), v.binSize()); }
};

CPSFIT_END_NAMESPACE
#endif
