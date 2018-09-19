#ifndef _PIPI_TO_SIGMA_RESAMPLED_CORRELATOR_H_
#define _PIPI_TO_SIGMA_RESAMPLED_CORRELATOR_H_

#include<config.h>
#include<utils/macros.h>

#include "data_containers.h"

CPSFIT_START_NAMESPACE

typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunction;

//Assumed to be a resampled distribution type
template<typename DistributionType>
correlationFunction<double, DistributionType> computePiPiToSigmaVacSub(const sigmaSelfContractionBase<DistributionType> &sigma_self, 
								       const bubbleDataBase<DistributionType> &pipi_self){

  int Lt = sigma_self.getLt(); assert(pipi_self.getLt() == Lt);
  int nsample = sigma_self(0).size(); assert(pipi_self(0).size() == nsample);

  correlationFunction<double, DistributionType> out(Lt, [&](const int t){ return typename correlationFunction<double, DistributionType>::ElementType(t, DistributionType(nsample,0.)); } );

  double coeff = -sqrt(6.)/2./double(Lt);

  for(int t0=0; t0<Lt; t0++){
    for(int tsep=0;tsep<Lt; tsep++){
      int t1 = (t0 + tsep) % Lt;
      
      out.value(tsep) = out.value(tsep) + coeff*pipi_self(t0)*sigma_self(t1);
    }
  }
  return out;
}

//Compute double-jackknife vacuum-subtraction from raw data
doubleJackCorrelationFunction computePiPiToSigmaVacSub(const sigmaSelfContraction &sigma_self, const bubbleData &pipi_self, const int bin_size){
  int Lt = sigma_self.getLt(); assert(pipi_self.getLt() == Lt);
  int nsample_raw = sigma_self(0).size(); assert(pipi_self(0).size() == nsample_raw);
  int nsample_binned = nsample_raw/bin_size;
  sigmaSelfContractionDoubleJack sigma_self_dj(Lt, nsample_binned);
  bubbleDataDoubleJack pipi_self_dj(Source, Lt, pipi_self.getTsepPiPi(), nsample_binned);
  for(int t=0;t<Lt;t++){
    sigma_self_dj(t).resample(sigma_self(t).bin(bin_size));
    pipi_self_dj(t).resample(pipi_self(t).bin(bin_size));
  }
  return computePiPiToSigmaVacSub(sigma_self_dj, pipi_self_dj);
}


template<typename DistributionType>
correlationFunction<double, DistributionType> foldPiPiToSigma(const correlationFunction<double, DistributionType> &data, const int Lt, const int tsep_pipi){
  correlationFunction<double, DistributionType> out(data);
  for(int t=0;t<Lt;t++){
    int tp = (Lt - tsep_pipi - t    + Lt) % Lt;
    out.value(t) = 0.5*( data.value(t) + data.value(tp) );
  }
  return out;
}

CPSFIT_END_NAMESPACE

#endif
