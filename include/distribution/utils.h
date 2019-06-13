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

CPSFIT_END_NAMESPACE
#endif
