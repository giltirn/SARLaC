#ifndef _CPSFIT_DISTRIBUTION_UTILS
#define _CPSFIT_DISTRIBUTION_UTILS

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>

#include<distribution/block_double_jackknife.h>
#include<distribution/bootstrap.h>
#include<distribution/boot_jackknife.h>


CPSFIT_START_NAMESPACE

template<typename ResampledDistributionType, typename RawDistributionType>
struct _bin_resample_t{ 
  inline static ResampledDistributionType doit(const RawDistributionType &raw, const int bin_size){ return RawDistributionType(raw.bin(bin_size)); } 
};
template<typename T, template<typename> class V, typename RawDistributionType>
struct _bin_resample_t<blockDoubleJackknifeDistribution<T,V>,  RawDistributionType>{ 
  inline static blockDoubleJackknifeDistribution<T,V> doit(const RawDistributionType &raw, const int bin_size){ return blockDoubleJackknifeDistribution<T,V>(raw, bin_size); }
};

template<typename ResampledDistributionType, typename RawDistributionType, 
	 typename std::enable_if<hasSampleMethod<ResampledDistributionType>::value && hasSampleMethod<RawDistributionType>::value, int>::type = 0
	>
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


template<typename ResampledDistributionType, typename RawDistributionType, typename binResampler, 
	 typename std::enable_if<hasSampleMethod<ResampledDistributionType>::value && hasSampleMethod<RawDistributionType>::value && !std::is_same<binResampler,int>::value,int>::type = 0>
inline ResampledDistributionType binResample(const RawDistributionType &raw, const binResampler &resampler){
  ResampledDistributionType out; resampler.binResample(out, raw); return out;
}



//Autocorrelation and integrated autocorrelation as defined in https://arxiv.org/pdf/1208.4412.pdf  page 9
double autocorrelation(const int delta, const rawDataDistribution<double> &data){
  double var = data.variance();
  double mean = data.mean();
  
  double c_delta = 0.;
  int navg = data.size()-delta;
  for(int t=0;t<data.size()-delta;t++){
    c_delta += ( data.sample(t) - mean ) * ( data.sample(t + delta) - mean ) / var;
  }
  c_delta /= navg;
  return c_delta;
}

//tau_int as a function of the cut on the separation.
double integratedAutocorrelation(const int delta_cut, const rawDataDistribution<double> &data){
  double tau_int = 0.5;
  for(int delta = 1; delta <= delta_cut; delta++)
    tau_int += autocorrelation(delta, data);
  return tau_int;
}

//Same as above but return result for every delta_cut from 0 .. delta_cut_max
std::vector<double> integratedAutocorrelationMulti(const int delta_cut_max, const rawDataDistribution<double> &data){
  double tau_int = 0.5;
    
  std::vector<double> out(1, tau_int);

  for(int delta = 1; delta <= delta_cut_max; delta++){
    tau_int += autocorrelation(delta,data);
    out.push_back(tau_int);
  }

  return out;
}

//Binned variant described in https://arxiv.org/pdf/1411.7017.pdf  page 16, bottom. Provides error bars
//For a consistent resampling we need to have the number of binned samples equal for all delta. 
//We therefore must provide a delta_max, from which we can generate nsample - delta_max  values of  (Y_i - Ybar)( Y_{i+delta} - Ybar )
//Resample table 'rtable' should be of size  nboots * ( (nsample-delta_max)/ bin_size )
bootstrapDistribution<double> autocorrelation(const int delta, const int bin_size, const int delta_max, const rawDataDistribution<double> &data, const std::vector<std::vector<int> > &rtable){
  double mean = data.mean();
  
  int nvals = data.size() - delta_max;

  rawDataDistribution<double> vals(nvals, [&](const int s){  return ( data.sample(s) - mean ) * ( data.sample(s + delta) - mean );  });
  rawDataDistribution<double> vals_0(nvals, [&](const int s){  return pow( data.sample(s) - mean,2); });

  vals = vals.bin(bin_size, true);
  vals_0 = vals_0.bin(bin_size, true);

  bootstrapDistribution<double> vals_b(vals, rtable);
  bootstrapDistribution<double> vals0_b(vals_0, rtable);

  return bootstrapDistribution<double>(vals_b/vals0_b);
}

//tau_int as a function of the cut on the separation.
//delta_max>=delta_cut
//Resample table 'rtable' should be of size  nboots * ( (nsample-delta_max)/ bin_size )
bootstrapDistribution<double> integratedAutocorrelation(const int delta_cut, const int bin_size, const int delta_max, const rawDataDistribution<double> &data, const std::vector<std::vector<int> > &rtable){
  int nboot = rtable.size();
  bootstrapDistribution<double> tau_int(0.5, bootstrapInitType(nboot));
  for(int delta = 1; delta <= delta_cut; delta++)
    tau_int = tau_int + autocorrelation(delta, bin_size, delta_max, data, rtable);
  return tau_int;
}

//delta_max>=delta_cut_max
std::vector<bootstrapDistribution<double> > integratedAutocorrelationMulti(const int delta_cut_max, const int bin_size, const int delta_max, const rawDataDistribution<double> &data, const std::vector<std::vector<int> > &rtable){
  int nboot = rtable.size();

  bootstrapDistribution<double> tau_int(0.5, bootstrapInitType(nboot));
  std::vector<bootstrapDistribution<double> > out(1, tau_int);

  for(int delta = 1; delta <= delta_cut_max; delta++){
    tau_int = tau_int + autocorrelation(delta, bin_size, delta_max, data, rtable);
    out.push_back(tau_int);
  }
  return out;
}





//In this version we resample the raw data prior to computing   ( v[s] - mu ) ( v[s+delta] - mu )
//Use block resample table to deal with the autocorrelations. This version may require larger block size but is free from abiguities surrounding inconsistency of the resample table
bootstrapDistribution<double> autocorrelationRawResample(const int delta, const rawDataDistribution<double> &data, const std::vector<std::vector<int> > &rtable){
  int nboot = rtable.size();

  bootstrapInitType init(nboot);
  bootstrapDistribution<double> out(init);
  out.best() = autocorrelation(delta, data);

#pragma omp parallel for
  for(int b=0;b<nboot;b++){
    int nsample = rtable[b].size();
    rawDataDistribution<double> data_r(nsample);
    for(int i=0;i<nsample;i++)
      data_r.sample(i) = data.sample( rtable[b][i] );
    
    out.sample(b) = autocorrelation(delta, data_r);
  }
  return out;
}

//tau_int as a function of the cut on the separation.
bootstrapDistribution<double> integratedAutocorrelationRawResample(const int delta_cut, const rawDataDistribution<double> &data, const std::vector<std::vector<int> > &rtable){
  int nboot = rtable.size();

  bootstrapInitType init(nboot);
  bootstrapDistribution<double> out(init);
  out.best() = integratedAutocorrelation(delta_cut, data);

#pragma omp parallel for
  for(int b=0;b<nboot;b++){
    int nsample = rtable[b].size();
    rawDataDistribution<double> data_r(nsample);
    for(int i=0;i<nsample;i++)
      data_r.sample(i) = data.sample( rtable[b][i] );
    
    out.sample(b) = integratedAutocorrelation(delta_cut, data_r);
  }
  return out;
}

std::vector<bootstrapDistribution<double> > integratedAutocorrelationMultiRawResample(const int delta_cut_max, const rawDataDistribution<double> &data, const std::vector<std::vector<int> > &rtable){
  int nboot = rtable.size();

  bootstrapInitType init(nboot);
  std::vector<bootstrapDistribution<double> > out(delta_cut_max+1,  bootstrapDistribution<double>(init));
  std::vector<double> tau_multi_best = integratedAutocorrelationMulti(delta_cut_max, data);
  for(int i=0;i<=delta_cut_max;i++) out[i].best() = tau_multi_best[i];

#pragma omp parallel for
  for(int b=0;b<nboot;b++){
    int nsample = rtable[b].size();
    rawDataDistribution<double> data_r(nsample);
    for(int i=0;i<nsample;i++)
      data_r.sample(i) = data.sample( rtable[b][i] );
    
    std::vector<double> tau_multi_r = integratedAutocorrelationMulti(delta_cut_max, data_r);
    for(int i=0;i<=delta_cut_max;i++) out[i].sample(b) = tau_multi_r[i];
  }
  return out;
}







CPSFIT_END_NAMESPACE
#endif
