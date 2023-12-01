#ifndef _FIT_SIMPLE_BOOTSTRAP_FIT_H_
#define _FIT_SIMPLE_BOOTSTRAP_FIT_H_

#include <fit_simple/fitfunc_manager.h>

typedef correlationFunction<double,double> correlationFunctionD;

template<typename ArgsType, typename CMDlineType>
void fit(bootstrapDistribution<parameterVectorD> &params,
	 bootstrapDistributionD &chisq,
	 int &dof,
	 const bootstrapCorrelationFunctionD &data_b,
	 const bootJackknifeCorrelationFunctionD &data_bj,
	 const ArgsType &args, const CMDlineType &cmdline,
	 bool do_write_and_plot = true, const std::string &filename_append = ""){
  
  //Get the data in the fit range
  bootstrapCorrelationFunctionD data_b_inrange = getDataInRange(data_b, args.t_min, args.t_max);
  bootJackknifeCorrelationFunctionD data_bj_inrange = getDataInRange(data_bj, args.t_min, args.t_max);

  std::cout << "All data\n";
  for(int i=0;i<data_b.size();i++){
    std::cout << (int)data_b.coord(i) << " " << data_b.value(i) << std::endl;
  }

  std::cout << "Data in range\n";
  for(int i=0;i<data_b_inrange.size();i++){
    std::cout << (int)data_b_inrange.coord(i) << " " << data_b_inrange.value(i) << std::endl;
  }

  //Get the fit function manager
  FitFuncManagerBase::Options opt;
  opt.load_guess = cmdline.load_guess;
  opt.guess_file = cmdline.guess_file;

  std::unique_ptr< FitFuncManagerBase > fitfunc_manager = getFitFuncManager(args.fitfunc, args.Lt, args.t_min, args.t_max, opt);

  //Set up the minimizer
  MarquardtLevenbergParameters<double> minparams;
  minparams.verbose = true;

  if(cmdline.load_mlparams){
    parse(minparams, cmdline.mlparams_file);
    std::cout << "Loaded minimizer params: " << minparams << std::endl;
  }
  
  simpleFitWrapper<bootstrapDistributionD> fitter(*fitfunc_manager->getFitFunc(), MinimizerType::MarquardtLevenberg, minparams);

  //Generate the covariance matrix
  switch(args.covariance_strategy){
  case CovarianceStrategy::FrozenCorrelated:
    fitter.generateCovarianceMatrix(data_b_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::Correlated:
    fitter.generateCovarianceMatrix(data_bj_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::Uncorrelated:
    fitter.generateCovarianceMatrix(data_bj_inrange, CostType::Uncorrelated);
    break;
  default:
    assert(0);
  }

  if(cmdline.save_covariance_matrix) fitter.writeCovarianceMatrixHDF5(cmdline.save_covariance_matrix_file);

  //Do the fit
  parameterVectorD guess = fitfunc_manager->getGuess();

  auto binit = data_b.value(0).getInitializer();
  params = bootstrapDistribution<parameterVectorD>(guess, binit);
  chisq = bootstrapDistributionD(binit);
    
  bootstrapDistributionD chisq_per_dof(binit);
  
  fitter.fit(params, chisq, chisq_per_dof, dof, data_b_inrange);

  if(do_write_and_plot){
    int nsample = data_bj.value(0).best().size();
    
    bootstrapDistributionD pvalue_chisq(binit), pvalue_Tsq(binit);
    for(int s=0;s<iterate<bootstrapDistributionD>::size(chisq);s++){
      iterate<bootstrapDistributionD>::at(s, pvalue_chisq) = chiSquareDistribution::pvalue(dof, iterate<bootstrapDistributionD>::at(s, chisq));
      iterate<bootstrapDistributionD>::at(s, pvalue_Tsq) = TsquareDistribution::pvalue(iterate<bootstrapDistributionD>::at(s, chisq), dof, nsample-1);
    }
    
    std::cout << "Params: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
    std::cout << "Dof: " << dof << std::endl;
    std::cout << "P-value(chi^2): " << pvalue_chisq << std::endl;
    std::cout << "P-value(T^2): " << pvalue_Tsq << std::endl;
    
#ifdef HAVE_HDF5
    writeParamsStandard(chisq, stringize("chisq%s.hdf5", filename_append.c_str()));
    writeParamsStandard(chisq_per_dof, stringize("chisq_per_dof%s.hdf5", filename_append.c_str()));
    writeParamsStandard(params, stringize("params%s.hdf5", filename_append.c_str())); 
    writeParamsStandard(pvalue_chisq, stringize("pvalue_chisq%s.hdf5", filename_append.c_str()));
    writeParamsStandard(pvalue_Tsq, stringize("pvalue_Tsq%s.hdf5", filename_append.c_str()));
#endif

    fitfunc_manager->plot(data_b, params);
  }
}




template<typename ArgsType, typename CMDlineType>
void fitCentral(parameterVectorD &params,
		double &chisq,
		int &dof,
		const correlationFunctionD &data_cen,
		const jackknifeCorrelationFunctionD &data_j,
		const ArgsType &args, const CMDlineType &cmdline, MarquardtLevenbergParameters<double> const *minparams_in = NULL){

  //Get the data in the fit range
  correlationFunctionD data_cen_inrange = getDataInRange(data_cen, args.t_min, args.t_max);
  jackknifeCorrelationFunctionD data_j_inrange = getDataInRange(data_j, args.t_min, args.t_max);

  //Get the fit function manager
  FitFuncManagerBase::Options opt;
  opt.load_guess = cmdline.load_guess;
  opt.guess_file = cmdline.guess_file;

  std::unique_ptr< FitFuncManagerBase > fitfunc_manager = getFitFuncManager(args.fitfunc, args.Lt, args.t_min, args.t_max, opt);

  //Set up the minimizer
  MarquardtLevenbergParameters<double> minparams;
  if(minparams_in != NULL) minparams = *minparams_in;
  
  simpleSingleFitWrapper fitter(*fitfunc_manager->getFitFunc(), MinimizerType::MarquardtLevenberg, minparams);

  //Generate the covariance matrix
  switch(args.covariance_strategy){
  case CovarianceStrategy::Correlated:
  case CovarianceStrategy::FrozenCorrelated:
    fitter.generateCovarianceMatrix(data_j_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::Uncorrelated:
    fitter.generateCovarianceMatrix(data_j_inrange, CostType::Uncorrelated);
    break;
  default:
    assert(0);
  }

  //Do the fit
  double chisq_per_dof;
  fitter.fit(params, chisq, chisq_per_dof, dof, data_cen_inrange);
  
  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Dof: " << dof << std::endl;
}



#endif
