#ifndef _PIPI_SIGMA_SIM_FIT_RESAMPLED_CORRELATOR_H_
#define _PIPI_SIGMA_SIM_FIT_RESAMPLED_CORRELATOR_H_


#include<config.h>
#include<utils/macros.h>

#include<pipi_common/pipi_common.h>

CPSFIT_START_NAMESPACE

doubleJackCorrelationFunction computePiPi2ptVacSub(const bubbleDataAllMomenta &raw, const int bin_size, const int tsep_pipi, 
						   const std::vector<threeMomentum> &pion_mom,
						   const PiPiProjector proj_src = PiPiProjector::A1, const PiPiProjector proj_snk = PiPiProjector::A1){
  PiPiCorrelatorBasicSelector corr_select(proj_src, proj_snk,PiPiMomAllowed::All,{0,0,0});

  bubbleDataDoubleJackAllMomenta dj_bubble_data = binDoubleJackknifeResampleBubble(raw, bin_size);
  doubleJackCorrelationFunction A2_realavg_V_dj = computePiPi2ptFigureVprojectSourceAvg(dj_bubble_data,tsep_pipi,corr_select,pion_mom);
  return doubleJackCorrelationFunction(3. * A2_realavg_V_dj);
}


doubleJackCorrelationFunction binDoubleJackResample(const rawCorrelationFunction &raw, const int bin_size){
  return doubleJackCorrelationFunction(raw.size(),
				       [&](const int t){ 
					 return typename doubleJackCorrelationFunction::ElementType(t, doubleJackknifeDistributionD(raw.value(t).bin(bin_size)));
				       }
				       );
}

CPSFIT_END_NAMESPACE

#endif
