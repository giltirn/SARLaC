#pragma once

#include "../figure_data_container.h"
#include "../sigma_bubble_data_container.h"

SARLAC_START_NAMESPACE

//Construct disconnected part from  Re(  sigma_bubble * sigma_bubble ) as we did in the parallel calculation
//sigma_self_data_Z should be pre-projected
void reconstructSigma2ptDisconnected(figureData &disconn, const sigmaSelfContractionZ &sigma_self_data_Z){
  const int Lt = sigma_self_data_Z.getLt();
  const int nsample = sigma_self_data_Z(0).size();
  
  disconn.setup(Lt, rawDataDistributionD(nsample));
  for(int tsrc=0;tsrc<Lt;tsrc++){
    for(int tsep=0;tsep<Lt;tsep++){
      int tsnk = ( tsrc + tsep ) % Lt;
      rawDataDistribution<std::complex<double> > valz = 0.5 * sigma_self_data_Z(tsrc) * sigma_self_data_Z(tsnk);
      for(int s=0;s<nsample;s++) disconn(tsrc,tsep).sample(s) = valz.sample(s).real();
    }
  }
}


//Compute the disconnected part from Re( sigma_bubble ) * Re (sigma_bubble)   - this has a slightly better statistical error than the above
//sigma_self_data should be pre-projected
void reconstructSigma2ptDisconnected(figureData &disconn, const sigmaSelfContraction &sigma_self_data){
  const int Lt = sigma_self_data.getLt();
  const int nsample = sigma_self_data(0).size();

  disconn.setup(Lt);
  for(int tsrc=0;tsrc<Lt;tsrc++){
    for(int tsep=0;tsep<Lt;tsep++){
      int tsnk = ( tsrc + tsep ) % Lt;
      disconn(tsrc,tsep) =  0.5 * sigma_self_data(tsrc) * sigma_self_data(tsnk);
    }
  }
}

void reconstructSigma2ptConnected(figureData &conn, const figureData &full, const figureData &disconn){
  const int Lt = full.getLt();
  assert(disconn.getLt() == Lt);

  const int nsample = full(0,0).size();
  assert( disconn(0,0).size() == nsample);

  conn.setup(Lt, rawDataDistributionD(nsample,0.));
  for(int tsrc=0; tsrc< Lt; tsrc++)
    for(int tsep=0;tsep<Lt;tsep++)
      conn(tsrc, tsep) = full(tsrc, tsep) - disconn(tsrc, tsep);
}



SARLAC_END_NAMESPACE
