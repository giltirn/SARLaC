#ifndef _SIGMA_RESAMPLED_CORRELATOR_H_
#define _SIGMA_RESAMPLED_CORRELATOR_H_

typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunction;


doubleJackCorrelationFunction computeSigmaVacSub(const sigmaSelfContraction &raw, const int bin_size){
  int Lt = raw.getLt();
  NumericVector<doubleJackknifeDistributionD> self_dj(Lt, [&](const int t){ return doubleJackknifeDistributionD(raw(t).bin(bin_size)); });

  int nsample = self_dj(0).size();

  doubleJackCorrelationFunction out(Lt, [&](const int t){ return doubleJackCorrelationFunction::ElementType(t, doubleJackknifeDistributionD(nsample,0.)); } );

  for(int t0=0; t0<Lt; t0++){
    for(int tsep=0;tsep<Lt; tsep++){
      int t1 = (t0 + tsep) % Lt;
      
      out.value(tsep) = out.value(tsep) + 0.5 * self_dj(t0)*self_dj(t1) / double(Lt);
    }
  }
  return out;
}

#endif
