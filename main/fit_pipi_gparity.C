#include <fstream>
#include <algorithm>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common_defs.h>
#include <sstream>

#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/correlationfunction.h>
#include <fit_pipi_gparity/read_data.h>
#include <fit_pipi_gparity/main.h>

int main(int argc, char* argv[]){
  const std::string data_dir = "/home/ckelly/CPS/build/CPSfit/pipi_data";
  const int tsep_pipi = 4;
  const int Lt = 64;
  const int traj_start = 984;
  const int traj_inc = 4;
  const int traj_lessthan = 996;
  const int nsample = (traj_lessthan - traj_start)/traj_inc;
  
  figureDataAllMomenta raw_data;
  readFigure(raw_data, 'C', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  readFigure(raw_data, 'D', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  readFigure(raw_data, 'R', data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  
  bubbleDataAllMomenta raw_bubble_data;
  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);
  computeV(raw_data, raw_bubble_data, tsep_pipi);

  figureData A2_C = projectA2('C', raw_data);
  figureData A2_D = projectA2('D', raw_data);
  figureData A2_R = projectA2('R', raw_data);
  figureData A2_V = projectA2('V', raw_data);

  typedef correlationFunction<distributionD> rawCorrelationFunction;
  
  rawCorrelationFunction A2_realavg_C = realSourceAverage(A2_C);
  rawCorrelationFunction A2_realavg_D = realSourceAverage(A2_D);
  rawCorrelationFunction A2_realavg_R = realSourceAverage(A2_R);
  rawCorrelationFunction A2_realavg_V = realSourceAverage(A2_V);

  rawCorrelationFunction pipi_raw = 2*A2_realavg_D + A2_realavg_C - 6*A2_realavg_R + 3*A2_realavg_V;

  typedef correlationFunction<doubleJackknifeDistributionD> doubleJackCorrelationFunction;
  
  doubleJackCorrelationFunction pipi_dj(Lt,
					[&pipi_raw,nsample](const int t)
					{
					  typename doubleJackCorrelationFunction::ElementType out(t, doubleJackknifeDistributionD(nsample));
					  out.second.resample(pipi_raw.value(t));
					  return out;
					}
					);
  bubbleDataDoubleJackAllMomenta dj_bubble_data = doubleJackknifeResampleBubble(raw_bubble_data);
  figureDataDoubleJackAllMomenta dj_data;  
  computeV(dj_data, dj_bubble_data, tsep_pipi);
  figureDataDoubleJack A2_V_dj = projectA2('V', dj_data);
  doubleJackCorrelationFunction A2_realavg_V_dj = realSourceAverage(A2_V_dj);

  pipi_dj = fold(pipi_dj, tsep_pipi);
  A2_realavg_V_dj = fold(A2_realavg_V_dj, tsep_pipi);

  //basicPrint<> printer;
  //printer << "pipi_dj: " << pipi_dj << std::endl;
  //printer << "V_dj: " << A2_realavg_V_dj << std::endl;
  
  doubleJackCorrelationFunction pipi_dj_vacsubbed = pipi_dj - 3*A2_realavg_V_dj;

  typedef correlationFunction<jackknifeDistributionD> jackknifeCorrelationFunction;

  jackknifeCorrelationFunction pipi_j_vacsubbed(Lt,
						[&pipi_dj_vacsubbed,nsample](const int t)
						{
						  typename jackknifeCorrelationFunction::ElementType out(t, jackknifeDistributionD(nsample));
						  out.second = pipi_dj_vacsubbed.value(t).toJackknife();
						  return out;
						}
						); 
  return 0;
}
