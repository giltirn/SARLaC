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
#include <fit_wrapper.h>
#include <sstream>
#include <boost/timer/timer.hpp>

#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>
#include <fit_pipi_gparity/fitfunc.h>
#include <fit_pipi_gparity/cmdline.h>
#include <fit_pipi_gparity/main.h>

int main(const int argc, const char* argv[]){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }
  
  CMDline cmdline(argc,argv,2);

  typedef FitCoshPlusConstant FitFunc;
  //typedef FitCoshPlusConstantDoubleExp FitFunc;

  
  FitFunc::Params guess;
  if(cmdline.load_guess){
    parse(guess, cmdline.guess_file);
    std::cout << "Loaded guess: " << guess << std::endl;
  }

  const int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  figureDataAllMomenta raw_data;
  bubbleDataAllMomenta raw_bubble_data;

  if(cmdline.load_data_checkpoint){
    loadCheckpoint<boost::archive::binary_iarchive>(raw_data, raw_bubble_data, cmdline.load_data_checkpoint_file);
  }else if(cmdline.load_text_data_checkpoint){
    loadCheckpoint<boost::archive::text_iarchive>(raw_data, raw_bubble_data, cmdline.load_text_data_checkpoint_file);    
  }else{
    readFigure(raw_data, 'C', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    readFigure(raw_data, 'D', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    readFigure(raw_data, 'R', args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    readBubble(raw_bubble_data, args.data_dir, args.tsep_pipi, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
  }

  if(cmdline.save_data_checkpoint){
    saveCheckpoint<boost::archive::binary_oarchive>(raw_data, raw_bubble_data, cmdline.save_data_checkpoint_file);
  }
  if(cmdline.save_text_data_checkpoint){
    saveCheckpoint<boost::archive::text_oarchive>(raw_data, raw_bubble_data, cmdline.save_text_data_checkpoint_file);
  }
  
  
  //Some of Daiqian's old data was measured on every source timeslice while the majority was measured every 8. To fix this discrepancy we explicitly zero the abnormal data  
  zeroUnmeasuredSourceTimeslices(raw_data, 'C', args.tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'D', args.tstep_pipi);
  zeroUnmeasuredSourceTimeslices(raw_data, 'R', args.tstep_pipi);

  computeV(raw_data, raw_bubble_data, args.tsep_pipi);

  figureData A2_C = projectA2('C', raw_data);
  figureData A2_D = projectA2('D', raw_data);
  figureData A2_R = projectA2('R', raw_data);
  figureData A2_V = projectA2('V', raw_data);

  typedef correlationFunction<double,rawDataDistributionD> rawCorrelationFunction;
  
  rawCorrelationFunction A2_realavg_C = sourceAverage(A2_C);
  rawCorrelationFunction A2_realavg_D = sourceAverage(A2_D);
  rawCorrelationFunction A2_realavg_R = sourceAverage(A2_R);
  rawCorrelationFunction A2_realavg_V = sourceAverage(A2_V);

  outputRawData("raw_data_Cpart.dat", A2_realavg_C, 1.);
  outputRawData("raw_data_Dpart.dat", A2_realavg_D, 2.);
  outputRawData("raw_data_Rpart.dat", A2_realavg_R, -6.);
  outputRawData("raw_data_Vpart.dat", A2_realavg_V, 3.);
  
  rawCorrelationFunction pipi_raw = 2*A2_realavg_D + A2_realavg_C - 6*A2_realavg_R + 3*A2_realavg_V;

  std::cout << "Raw data:\n" << pipi_raw << std::endl;

  typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
  
  doubleJackCorrelationFunction pipi_dj(args.Lt,
					[&pipi_raw,nsample](const int t)
					{
					  typename doubleJackCorrelationFunction::ElementType out(t, doubleJackknifeDistributionD(nsample));
					  out.second.resample(pipi_raw.value(t));
					  return out;
					}
					);
  bubbleDataDoubleJackAllMomenta dj_bubble_data = doubleJackknifeResampleBubble(raw_bubble_data);
  doubleJackCorrelationFunction A2_realavg_V_dj = computeVprojectA2sourceAvg(dj_bubble_data,args.tsep_pipi);     //sourceAverage(A2_V_dj);
  
  pipi_dj = fold(pipi_dj, args.tsep_pipi);
  A2_realavg_V_dj = fold(A2_realavg_V_dj, args.tsep_pipi);

  doubleJackCorrelationFunction pipi_dj_vacsubbed = pipi_dj - 3*A2_realavg_V_dj;

  filterXrange<double> trange(args.t_min,args.t_max);
  
  filteredDataSeries<doubleJackCorrelationFunction> pipi_dj_vacsubbed_inrange(pipi_dj_vacsubbed, trange);

  const int ndata_fit = pipi_dj_vacsubbed_inrange.size();
  
  typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunction;

  jackknifeCorrelationFunction pipi_j_vacsubbed_inrange(ndata_fit,
							[&pipi_dj_vacsubbed_inrange](const int i)
							{
							  return typename jackknifeCorrelationFunction::ElementType(double(pipi_dj_vacsubbed_inrange.coord(i)), pipi_dj_vacsubbed_inrange.value(i).toJackknife());
							}
							);
  
  NumericSquareMatrix<jackknifeDistributionD> cov(ndata_fit);
  NumericVector<jackknifeDistributionD> sigma(ndata_fit);
  for(int i=0;i<ndata_fit;i++){
    cov(i,i) = doubleJackknifeDistributionD::covariance(pipi_dj_vacsubbed_inrange.value(i),  pipi_dj_vacsubbed_inrange.value(i));
    sigma(i) = sqrt(cov(i,i));
    
    for(int j=i+1;j<ndata_fit;j++)
      cov(i,j) = cov(j,i) = doubleJackknifeDistributionD::covariance(pipi_dj_vacsubbed_inrange.value(i),  pipi_dj_vacsubbed_inrange.value(j));
  }


  NumericSquareMatrix<jackknifeDistributionD> corr(ndata_fit);

  for(int i=0;i<ndata_fit;i++){
    corr(i,i) = jackknifeDistributionD(nsample,1.);
    
    for(int j=i+1;j<ndata_fit;j++)
      corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
  }

  NumericSquareMatrix<jackknifeDistributionD> inv_corr(ndata_fit, jackknifeDistributionD(nsample));
  svd_inverse(inv_corr, corr);
  std::cout << "Correlation matrix:\n" << corr << std::endl;
  std::cout << "Inverse correlation matrix:\n" << inv_corr << std::endl;

  NumericSquareMatrix<jackknifeDistributionD> test = corr * inv_corr;
  for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistributionD(nsample,1.0);

  std::cout << "|CorrMat * CorrMat^{-1} - 1|^2 = " << mod2(test) << std::endl;

  FitFunc fitfunc(args.Lt, args.tsep_pipi, args.Ascale, args.Cscale);

  typedef typename composeFitPolicy<double, FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  //typedef typename composeFitPolicy<double, FitFunc, frozenFitFuncPolicy, uncorrelatedFitPolicy>::type FitPolicies;
  fitter<FitPolicies> fitter;
  fitter.importFitFunc(fitfunc);
  fitter.importCostFunctionParameters(inv_corr,sigma);
  //fitter.importCostFunctionParameters(sigma);  

  //jackknifeDistribution<FitFunc::Params> freeze(nsample, FitFunc::Params(0,0,0,0,0));
  //fitter.freeze({4}, freeze);
  
  jackknifeDistribution<FitFunc::Params> params(nsample, guess);
  jackknifeDistributionD chisq;
  jackknifeDistributionD chisq_per_dof;
  fitter.fit(params, chisq, chisq_per_dof, pipi_j_vacsubbed_inrange);

  distributionPrint<decltype(params)>::printer(new pipiParamsPrinter<FitFunc>);

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
#endif
  
  return 0;
}
