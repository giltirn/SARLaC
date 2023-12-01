#ifndef _JACKKNIFE_C_DIST_PLOT_H_
#define _JACKKNIFE_C_DIST_PLOT_H_

#include<config.h>
#include<distribution/jackknife_C/class.h>
#include<plot/plot/accessors_2d.h>

SARLAC_START_NAMESPACE

template<typename T, template<typename> class V>
class DistributionPlotAccessor<jackknifeCdistribution<T,V> >{
public:
  static inline double value(const jackknifeCdistribution<T,V> &d){ return d.best(); }
  static inline double errplus(const jackknifeCdistribution<T,V> &d){ return d.standardError(); }
  static inline double errminus(const jackknifeCdistribution<T,V> &d){ return d.standardError(); }  
};

SARLAC_END_NAMESPACE
#endif
