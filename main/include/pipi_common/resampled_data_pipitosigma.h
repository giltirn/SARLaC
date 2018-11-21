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

  correlationFunction<double, DistributionType> out(Lt);

  double coeff = -sqrt(6.)/2./double(Lt);

  for(int tsep=0;tsep<Lt; tsep++){
    out.coord(tsep) = tsep;
    for(int t0=0; t0<Lt; t0++){
      int t1 = (t0 + tsep) % Lt;
      DistributionType tmp = coeff*pipi_self(t0)*sigma_self(t1);
      out.value(tsep) = t0 == 0 ? tmp : out.value(tsep) + tmp;
    }
  }
  return out;
}

//Compute double-jackknife vacuum-subtraction from raw data
template<typename resampledCorrelationFunctionType>
resampledCorrelationFunctionType computePiPiToSigmaVacSub(const sigmaSelfContraction &sigma_self, const bubbleData &pipi_self, const int bin_size){
  int Lt = sigma_self.getLt(); assert(pipi_self.getLt() == Lt);
  int nsample_raw = sigma_self(0).size(); assert(pipi_self(0).size() == nsample_raw);

  typedef typename resampledCorrelationFunctionType::DataType DistributionType;
  typedef sigmaSelfContractionSelect<DistributionType> sigmaSelfContractionType;
  typedef bubbleDataSelect<DistributionType> pipiSelfContractionType;

  sigmaSelfContractionType sigma_self_r(Lt);
  pipiSelfContractionType pipi_self_r(Source, Lt, pipi_self.getTsepPiPi());
  for(int t=0;t<Lt;t++){
    sigma_self_r(t).resample(sigma_self(t).bin(bin_size));
    pipi_self_r(t).resample(pipi_self(t).bin(bin_size));
  }
  return computePiPiToSigmaVacSub(sigma_self_r, pipi_self_r);
}
template<typename resampledCorrelationFunctionType>
inline resampledCorrelationFunctionType computePiPiToSigmaVacSub(const sigmaSelfContraction &sigma_self, const bubbleDataAllMomenta &pipi_self, const PiPiProjectorBase &proj_pipi, const int bin_size){
  auto proj_bubble = projectSourcePiPiBubble(pipi_self, proj_pipi);
  return computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(sigma_self, proj_bubble , bin_size);
}
template<typename resampledCorrelationFunctionType>
inline resampledCorrelationFunctionType computePiPiToSigmaVacSub(const sigmaSelfContraction &sigma_self, const bubbleDataAllMomenta &pipi_self, const PiPiProjector proj_pipi, const int bin_size){
  std::unique_ptr<PiPiProjectorBase> proj(getProjector(proj_pipi,{0,0,0}));
  return computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(sigma_self, pipi_self, *proj, bin_size);
}


CPSFIT_END_NAMESPACE

#endif
