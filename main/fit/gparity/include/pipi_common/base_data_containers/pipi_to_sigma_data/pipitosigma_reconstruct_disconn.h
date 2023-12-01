#pragma once

#include<config.h>
#include<utils/macros.h>

#include "pipitosigma_data.h"

SARLAC_START_NAMESPACE

//Construct disconnected part from  Re(  pipi_buble * sigma_bubble ) as we did in the parallel calculation
//pipi_self_data_Z should be pre-projected
void reconstructPiPiToSigmaDisconnected(figureData &disconn, const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z){
  const int Lt = pipi_self_data_Z.getLt();
  assert(sigma_self_data_Z.getLt() == Lt);

  const int nsample = pipi_self_data_Z(0).size();
  assert( sigma_self_data_Z(0).size() == nsample);

  disconn.setup(Lt, rawDataDistributionD(nsample));
  for(int tsrc=0;tsrc<Lt;tsrc++){
    for(int tsep=0;tsep<Lt;tsep++){
      int tsnk = ( tsrc + tsep ) % Lt;
      //note the pipi bubble is computed online as 0.5 tr( mf(t) mf(t-tsep) ), and this is combined with the sigma bubble with a coeff sqrt(6)/2 in the parallel code to form the pipi->sigma disconnected part
      //However the correct formula for the pipi bubble is -0.5 tr( mf(t) mf(t-tsep) )   ; this is corrected for when the bubble is loaded in this analysis code. Hence the coeff here needs to be -sqrt(6)/2
      rawDataDistribution<std::complex<double> > valz = -sqrt(6.)/2 * pipi_self_data_Z(tsrc) * sigma_self_data_Z(tsnk); 
      for(int s=0;s<nsample;s++) disconn(tsrc,tsep).sample(s) = valz.sample(s).real();
    }
  }
}

//Compute the disconnected part from Re( pipi_bubble ) * Re (sigma_bubble)   - this has a slightly better statistical error than the above
//pipi_self_data should be pre-projected
void reconstructPiPiToSigmaDisconnected(figureData &disconn, const bubbleData &pipi_self_data, const sigmaSelfContraction &sigma_self_data){
  const int Lt = pipi_self_data.getLt();
  assert(sigma_self_data.getLt() == Lt);

  const int nsample = pipi_self_data(0).size();
  assert( sigma_self_data(0).size() == nsample);

  disconn.setup(Lt);
  for(int tsrc=0;tsrc<Lt;tsrc++){
    for(int tsep=0;tsep<Lt;tsep++){
      int tsnk = ( tsrc + tsep ) % Lt;
      disconn(tsrc,tsep) = -sqrt(6.)/2 * pipi_self_data(tsrc) * sigma_self_data(tsnk); 
    }
  }
}



//In the current runs we sum the connected and disconnected diagrams online. However for tstep > 1 this means the vacuum diagrams are also only sampled on a subset of timeslices
//This can be rectified by removing the disconnected component and performing the source timeslice average of the connected part, then afterwards constructing the disconnected part
//and source timeslice averaging that over all Lt
void reconstructPiPiToSigmaConnected(figureData &conn, const figureData &full, const figureData &disconn, const int tstep_src_full){
  const int Lt = full.getLt();
  assert(disconn.getLt() == Lt);

  const int nsample = full(0,0).size();
  assert( disconn(0,0).size() == nsample);

  conn.setup(Lt, rawDataDistributionD(nsample,0.));
  for(int tsrc=0; tsrc< Lt; tsrc += tstep_src_full)
    for(int tsep=0;tsep<Lt;tsep++)
      conn(tsrc, tsep) = full(tsrc, tsep) - disconn(tsrc, tsep);
}



SARLAC_END_NAMESPACE
