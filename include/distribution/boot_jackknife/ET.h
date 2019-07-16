#ifndef _BOOT_JACKKNIFE_DIST_ET_H_
#define _BOOT_JACKKNIFE_DIST_ET_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/boot_jackknife/class.h>

CPSFIT_START_NAMESPACE

template<typename A, template<typename> class V>
  struct getElem<bootJackknifeDistribution<A,V> >{
  static inline auto elem(const bootJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline auto elem(bootJackknifeDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return v.sample(i); }
  static inline bootJackknifeInitType common_properties(const bootJackknifeDistribution<A,V> &v){ return v.getInitializer(); }
};

CPSFIT_END_NAMESPACE
#endif
