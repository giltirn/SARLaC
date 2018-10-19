#ifndef _SIGMA_RESAMPLED_CORRELATOR_H_
#define _SIGMA_RESAMPLED_CORRELATOR_H_

#include<config.h>
#include<utils/macros.h>

#include "data_containers.h"

CPSFIT_START_NAMESPACE

typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunction;

//Assumed to be a resampled distribution type
template<typename DistributionType>
correlationFunction<double, DistributionType> computeSigmaVacSub(const sigmaSelfContractionBase<DistributionType> &bub){
  int nsample = bub(0).size();
  int Lt = bub.getLt();

  correlationFunction<double,DistributionType> out(Lt, [&](const int t){ return typename correlationFunction<double,DistributionType>::ElementType(t, DistributionType(nsample,0.)); } );

  for(int t0=0; t0<Lt; t0++){
    for(int tsep=0;tsep<Lt; tsep++){
      int t1 = (t0 + tsep) % Lt;
      
      out.value(tsep) = out.value(tsep) + 0.5 * bub(t0)*bub(t1) / double(Lt);
    }
  }
  return out;
}

template<typename resampledCorrelationFunctionType>
resampledCorrelationFunctionType computeSigmaVacSub(const sigmaSelfContraction &raw, const int bin_size){
  int Lt = raw.getLt();
  int nsample_raw = raw(0).size();
  int nsample_binned = nsample_raw/bin_size;
  
  typedef typename resampledCorrelationFunctionType::DataType DistributionType;
  typedef sigmaSelfContractionSelect<DistributionType> sigmaSelfContractionType;
  sigmaSelfContractionType sigma_self_r(Lt, nsample_binned);
  for(int t=0;t<Lt;t++)
    sigma_self_r(t).resample(raw(t).bin(bin_size));
  return computeSigmaVacSub(sigma_self_r);
}

CPSFIT_END_NAMESPACE

#endif
