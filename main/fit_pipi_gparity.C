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
#include <fit_pipi_gparity/fitfunc.h>
#include <fit_pipi_gparity/cmdline.h>
#include <fit_pipi_gparity/main.h>

int main(const int argc, const char* argv[]){
  CMDline cmdline(argc,argv,1);
  
  const std::string data_dir = "/home/ckelly/CPS/build/CPSfit/pipi_data";
  const int tsep_pipi = 4;
  const int Lt = 64;
  const int traj_start = 984;
  const int traj_inc = 4;
  const int traj_lessthan = 996;
  const int nsample = (traj_lessthan - traj_start)/traj_inc;
  const int t_min = 6;
  const int t_max = 25;
  const double Ascale = 1e13;
  const double Cscale = 1e13;

  typedef FitCoshPlusConstant FitFunc;
  FitFunc::Params guess;
  if(cmdline.load_guess){
    std::ifstream f(cmdline.guess_file.c_str());
    assert(f.good());
    f >> guess;
    assert(!f.bad() && !f.fail());
    f.close();

    std::cout << "Loaded guess: " << guess << std::endl;
  }else{
    guess.A = 1;
    guess.E = 0.3;
    guess.C = 0;
  }
  
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

  distributionPrint<doubleJackknifeDistributionD>::printer(new centralValueDistributionPrinter<doubleJackknifeDistributionD>);
  std::cout << "pipi_dj: " << pipi_dj << std::endl;
  std::cout << "V_dj: " << A2_realavg_V_dj << std::endl;
  
  doubleJackCorrelationFunction pipi_dj_vacsubbed = pipi_dj - 3*A2_realavg_V_dj;

  filterXrange<double> trange(t_min,t_max);
  
  filteredDataSeries<doubleJackCorrelationFunction> pipi_dj_vacsubbed_inrange(pipi_dj_vacsubbed, trange);

  const int ndata_fit = pipi_dj_vacsubbed_inrange.size();
  
  typedef correlationFunction<jackknifeDistributionD> jackknifeCorrelationFunction;

  jackknifeCorrelationFunction pipi_j_vacsubbed_inrange(ndata_fit,
							[&pipi_dj_vacsubbed_inrange](const int i)
							{
							  return typename jackknifeCorrelationFunction::ElementType(double(pipi_dj_vacsubbed_inrange.coord(i)), pipi_dj_vacsubbed_inrange.value(i).toJackknife());
							}
							);
  
  NumericMatrix<jackknifeDistributionD> cov(ndata_fit);
  NumericVector<jackknifeDistributionD> sigma(ndata_fit);
  for(int i=0;i<ndata_fit;i++){
    cov(i,i) = doubleJackknifeDistributionD::covariance(pipi_dj_vacsubbed_inrange.value(i),  pipi_dj_vacsubbed_inrange.value(i));
    sigma(i) = sqrt(cov(i,i));
    
    for(int j=i+1;j<ndata_fit;j++)
      cov(i,j) = cov(j,i) = doubleJackknifeDistributionD::covariance(pipi_dj_vacsubbed_inrange.value(i),  pipi_dj_vacsubbed_inrange.value(j));
  }


  NumericMatrix<jackknifeDistributionD> corr(ndata_fit);

  for(int i=0;i<ndata_fit;i++){
    corr(i,i) = jackknifeDistributionD(nsample,1.);
    
    for(int j=i+1;j<ndata_fit;j++)
      corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
  }

  NumericMatrix<jackknifeDistributionD> inv_corr(ndata_fit, jackknifeDistributionD(nsample));
  svd_inverse(inv_corr, corr);

  FitFunc fitfunc(Lt, tsep_pipi, Ascale, Cscale);

  jackknifeDistribution<FitFunc::Params> params(nsample, guess);
  jackknifeDistributionD chisq(nsample);
  
  typedef sampleSeries<const decltype(pipi_j_vacsubbed_inrange)> sampleSeriesType;
  typedef NumericMatrixSampleView<const decltype(inv_corr)> sampleInvCorrType;
  typedef CorrelatedChisqCostFunction<FitFunc, sampleSeriesType, sampleInvCorrType> costFunctionType;
  typedef MarquardtLevenbergMinimizer<costFunctionType> minimizerType;
  typedef minimizerType::AlgorithmParameterType minimizerParamsType;

  minimizerParamsType min_params;
  
#pragma omp parallel for
  for(int s=0;s<nsample;s++){
    sampleSeriesType data_s(pipi_j_vacsubbed_inrange, s);
    sampleInvCorrType inv_corr_s(inv_corr, s);
  
    std::vector<double> sigma_s(ndata_fit);
    for(int d=0;d<ndata_fit;d++) sigma_s[d] = sigma[d].sample(s);
    
    costFunctionType cost_func(fitfunc, data_s, sigma_s, inv_corr_s);
    minimizerType minimizer(cost_func,min_params);
    assert(minimizer.hasConverged());
  }

  int dof = ndata_fit - fitfunc.Nparams();
  jackknifeDistributionD chisq_per_dof = chisq/double(dof);

  distributionPrint<decltype(params)>::printer(new pipiParamsPrinter);

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
  return 0;
}
