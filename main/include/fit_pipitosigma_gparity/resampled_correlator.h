#ifndef _PIPI_TO_SIGMA_RESAMPLED_CORRELATOR_H_
#define _PIPI_TO_SIGMA_RESAMPLED_CORRELATOR_H_

typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunction;


doubleJackCorrelationFunction computePiPiToSigmaVacSub(const sigmaSelfContraction &sigma_self, const bubbleData &pipi_self, const int bin_size){
  int Lt = sigma_self.getLt();
  NumericVector<doubleJackknifeDistributionD> sigma_self_dj(Lt, [&](const int t){ return doubleJackknifeDistributionD(sigma_self(t).bin(bin_size)); });
  NumericVector<doubleJackknifeDistributionD> pipi_self_dj(Lt, [&](const int t){ return doubleJackknifeDistributionD(pipi_self(t).bin(bin_size)); });

  int nsample = sigma_self_dj(0).size();

  doubleJackCorrelationFunction out(Lt, [&](const int t){ return doubleJackCorrelationFunction::ElementType(t, doubleJackknifeDistributionD(nsample,0.)); } );

  double coeff = -sqrt(6.)/2./double(Lt);

  for(int t0=0; t0<Lt; t0++){
    for(int tsep=0;tsep<Lt; tsep++){
      int t1 = (t0 + tsep) % Lt;
      
      out.value(tsep) = out.value(tsep) + coeff*pipi_self_dj(t0)*sigma_self_dj(t1);
    }
  }
  return out;
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

#endif
