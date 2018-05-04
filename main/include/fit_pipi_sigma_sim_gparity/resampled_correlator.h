#ifndef _PIPI_SIGMA_SIM_FIT_RESAMPLED_CORRELATOR_H_
#define _PIPI_SIGMA_SIM_FIT_RESAMPLED_CORRELATOR_H_

doubleJackCorrelationFunction computePiPi2ptVacSub(const bubbleDataAllMomenta &raw, const int bin_size, const int tsep_pipi, const std::vector<threeMomentum> &pion_mom){
  PiPiProjectA1 proj;
  PiPiMomAllowAll allow;
  bubbleDataDoubleJackAllMomenta dj_bubble_data = binDoubleJackknifeResampleBubble(raw, bin_size);
  doubleJackCorrelationFunction A2_realavg_V_dj = computeVprojectSourceAvg(dj_bubble_data,tsep_pipi,proj,proj,allow,pion_mom);
  return doubleJackCorrelationFunction(3. * A2_realavg_V_dj);
}


doubleJackCorrelationFunction binDoubleJackResample(const rawCorrelationFunction &raw, const int bin_size){
  return doubleJackCorrelationFunction(raw.size(),
				       [&](const int t){ 
					 return typename doubleJackCorrelationFunction::ElementType(t, doubleJackknifeDistributionD(raw.value(t).bin(bin_size)));
				       }
				       );
}



#endif
