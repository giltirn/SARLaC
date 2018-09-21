#ifndef _PIPI_SIGMA_SIM_FIT_RESAMPLED_CORRELATOR_H_
#define _PIPI_SIGMA_SIM_FIT_RESAMPLED_CORRELATOR_H_


#include<config.h>
#include<utils/macros.h>

#include<pipi_common/pipi_common.h>

CPSFIT_START_NAMESPACE

doubleJackCorrelationFunction computePiPi2ptVacSub(const bubbleDataAllMomenta &raw, const int bin_size, const int tsep_pipi, 
						   const std::vector<threeMomentum> &pion_mom,
						   const PiPiProjector proj_src_t = PiPiProjector::A1momSet111, const PiPiProjector proj_snk_t = PiPiProjector::A1momSet111){
  std::unique_ptr<PiPiProjectorBase> proj_src( getProjector(proj_src_t, {0,0,0}) );
  std::unique_ptr<PiPiProjectorBase> proj_snk( getProjector(proj_snk_t, {0,0,0}) );

  bubbleDataDoubleJackAllMomenta dj_bubble_data = binDoubleJackknifeResampleBubble(raw, bin_size);
  doubleJackCorrelationFunction A2_realavg_V_dj = computePiPi2ptFigureVprojectSourceAvg(dj_bubble_data,tsep_pipi,*proj_src,*proj_snk);
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
