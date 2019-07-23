#ifndef _CPSFIT_DISTRIBUTION_UTILS
#define _CPSFIT_DISTRIBUTION_UTILS

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>

#include<distribution/block_double_jackknife.h>

CPSFIT_START_NAMESPACE

template<typename ResampledDistributionType, typename RawDistributionType>
struct _bin_resample_t{ 
  inline static ResampledDistributionType doit(const RawDistributionType &raw, const int bin_size){ return RawDistributionType(raw.bin(bin_size)); } 
};
template<typename T, template<typename> class V, typename RawDistributionType>
struct _bin_resample_t<blockDoubleJackknifeDistribution<T,V>,  RawDistributionType>{ 
  inline static blockDoubleJackknifeDistribution<T,V> doit(const RawDistributionType &raw, const int bin_size){ return blockDoubleJackknifeDistribution<T,V>(raw, bin_size); }
};

template<typename ResampledDistributionType, typename RawDistributionType>
ResampledDistributionType binResample(const RawDistributionType &raw, const int bin_size){
  return _bin_resample_t<ResampledDistributionType, RawDistributionType>::doit(raw, bin_size);
}

//These resampler classes can be used to generalized the concepts of binning and resampling in code that computes resampled data from raw
struct basicBinResampler{
  int bin_size;
  basicBinResampler(int bin_size): bin_size(bin_size){}
  basicBinResampler(): bin_size(0){}

  template<template<typename,template<typename> class> class DistributionType, typename T, template<typename> class V>
  inline void binResample(DistributionType<T,V> &out, const rawDataDistribution<T,V> &in) const{ out.resample(in.bin(bin_size)); }
  
  template<typename T, template<typename> class V>
  inline void binResample(blockDoubleJackknifeDistribution<T,V> &out, const rawDataDistribution<T,V> &in) const{ out.resample(in, bin_size); }

  template<typename T, template<typename> class V>
  inline void binResample(T &out, const rawDataDistribution<T,V> &in) const{  out = in.bin(bin_size).mean(); }
};

//We don't bin, rather the resample table should use a block resampling strategy
struct bootstrapBlockResampler{
  const std::vector<std::vector<int> > &rtable;

  bootstrapBlockResampler(const std::vector<std::vector<int> > &rtable): rtable(rtable){}

  template<typename T, template<typename> class V>
  inline void binResample(bootstrapDistribution<T,V> &out, const rawDataDistribution<T,V> &in) const{ out.resample(in, rtable); }
  template<typename T, template<typename> class V>
  inline void binResample(bootJackknifeDistribution<T,V> &out, const rawDataDistribution<T,V> &in) const{ out.resample(in, rtable); }
};


template<typename ResampledDistributionType, typename RawDistributionType, typename binResampler>
inline ResampledDistributionType binResample(const RawDistributionType &raw, const binResampler &resampler){
  ResampledDistributionType out; resampler.binResample(out, raw); return out;
}

CPSFIT_END_NAMESPACE
#endif
