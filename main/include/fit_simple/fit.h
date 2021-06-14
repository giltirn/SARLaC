#ifndef _FIT_SIMPLE_FIT_H_
#define _FIT_SIMPLE_FIT_H_

#include "fitfunc_manager.h"

template<typename ArgsType, typename CMDlineType>
void fit(jackknifeDistribution<parameterVectorD> &params,
	 jackknifeDistributionD &chisq,
	 int &dof,
	 const jackknifeCorrelationFunctionD &data_j,
	 const doubleJackknifeCorrelationFunctionD &data_dj,
	 const blockDoubleJackknifeCorrelationFunctionD &data_bdj,
	 const ArgsType &args, const CMDlineType &cmdline){
  
  bool do_dj, do_bdj;
  getDJtypes(do_dj, do_bdj, args.covariance_strategy);

  //Get the data in the fit range
  jackknifeCorrelationFunctionD data_j_inrange = getDataInRange(data_j, args.t_min, args.t_max);
  doubleJackknifeCorrelationFunctionD data_dj_inrange;
  if(do_dj) data_dj_inrange = getDataInRange(data_dj, args.t_min, args.t_max);
  blockDoubleJackknifeCorrelationFunctionD data_bdj_inrange;
  if(do_bdj) data_bdj_inrange = getDataInRange(data_bdj, args.t_min, args.t_max);

  std::cout << "All data\n";
  for(int i=0;i<data_j.size();i++){
    std::cout << (int)data_j.coord(i) << " " << data_j.value(i) << std::endl;
  }

  std::cout << "Data in range\n";
  for(int i=0;i<data_j_inrange.size();i++){
    std::cout << (int)data_j_inrange.coord(i) << " " << data_j_inrange.value(i) << std::endl;
  }


  //Get the fit function manager
  FitFuncManagerBase::Options opt;
  opt.load_guess = cmdline.load_guess;
  opt.guess_file = cmdline.guess_file;

  std::unique_ptr< FitFuncManagerBase > fitfunc_manager = getFitFuncManager(args.fitfunc, args.Lt, args.t_min, args.t_max, opt);

  //Set up the minimizer
  MarquardtLevenbergParameters<double> minparams;
  if(cmdline.load_mlparams){
    parse(minparams, cmdline.mlparams_file);
    std::cout << "Loaded minimizer params: " << minparams << std::endl;
  }else{
    minparams.verbose = true;
  }
  
  simpleFitWrapper<jackknifeDistributionD> fitter(*fitfunc_manager->getFitFunc(), MinimizerType::MarquardtLevenberg, minparams);

  //Generate the covariance matrix
  switch(args.covariance_strategy){
  case CovarianceStrategy::CorrelatedBlockHybrid:
    fitter.generateCovarianceMatrix(data_dj_inrange, data_bdj_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::FrozenCorrelated:
    fitter.generateCovarianceMatrix(data_j_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::Correlated:
    fitter.generateCovarianceMatrix(data_dj_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::CorrelatedBlock:
    fitter.generateCovarianceMatrix(data_bdj_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::Uncorrelated:
    fitter.generateCovarianceMatrix(data_dj_inrange, CostType::Uncorrelated);
    break;
  default:
    assert(0);
  }

  //Do the fit
  parameterVectorD guess = fitfunc_manager->getGuess();

  const int nsample = data_j.value(0).size();
  params = jackknifeDistribution<parameterVectorD>(nsample, guess);
  chisq = jackknifeDistributionD(nsample);
    
  jackknifeDistributionD chisq_per_dof(nsample);
  
  fitter.fit(params, chisq, chisq_per_dof, dof, data_j_inrange);

  std::cout << "Completed fit, chisq: " << chisq << std::endl;

  gsl_error_handler_t * err = gsl_set_error_handler_off();
  jackknifeDistributionD pvalue_chisq(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });
  jackknifeDistributionD pvalue_Tsq(nsample, 0.);
  bool do_Tsq = nsample - dof >= 0;
  if(!do_Tsq) std::cout << "WARNING: Skipping T^2 p-value because N-dof=" << nsample-dof << "<0" << std::endl;

  if(do_Tsq) pvalue_Tsq = jackknifeDistributionD(nsample, [&](const int s){ return TsquareDistribution::pvalue(chisq.sample(s), dof, nsample-1); });

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "Dof: " << dof << std::endl;
  std::cout << "P-value(chi^2): " << pvalue_chisq << std::endl;
  if(do_Tsq) std::cout << "P-value(T^2): " << pvalue_Tsq << std::endl;
  
#ifdef HAVE_HDF5
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(params, "params.hdf5"); 
  writeParamsStandard(pvalue_chisq, "pvalue_chisq.hdf5");
  if(do_Tsq) writeParamsStandard(pvalue_Tsq, "pvalue_Tsq.hdf5");
#endif

  fitfunc_manager->plot(data_j, params);
}



template<typename ArgsType, typename CMDlineType>
void fitCentral(parameterVectorD &params,
		double &chisq,
		int &dof,
		const jackknifeCorrelationFunctionD &data_j,
		const jackknifeCorrelationFunctionD &data_j_unbinned,
		const ArgsType &args, const CMDlineType &cmdline, MarquardtLevenbergParameters<double> const *minparams_in = NULL){

  //For block and block-hybrid the unbinned jackknife is used, for the former to compute the covariance matrix 
  //and for the latter the correlation matrix (with sigma computed from the binned jackknife)
  bool do_j_b, do_j_ub;
  getJtypes(do_j_b, do_j_ub, args.covariance_strategy);

  //Get the data in the fit range
  jackknifeCorrelationFunctionD data_j_inrange, data_j_ub_inrange;
  if(do_j_b) data_j_inrange = getDataInRange(data_j, args.t_min, args.t_max);
  if(do_j_ub) data_j_ub_inrange = getDataInRange(data_j_unbinned, args.t_min, args.t_max);

  //Get the central values
  const jackknifeCorrelationFunctionD & data_j_inrange_use = do_j_b ?  data_j_inrange : data_j_ub_inrange;

  correlationFunction<double, double> data_cen_inrange(data_j_inrange_use.size());  
  for(int i=0;i<data_j_inrange_use.size();i++){
    data_cen_inrange.coord(i) = data_j_inrange_use.coord(i);
    data_cen_inrange.value(i) = data_j_inrange_use.value(i).mean();
  }

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
  case CovarianceStrategy::CorrelatedBlockHybrid:
    fitter.generateCovarianceMatrix(data_j_inrange, data_j_ub_inrange);
    break;
  case CovarianceStrategy::FrozenCorrelated:
    //Freezing the covariance matrix doesn't make any difference if only one sample!
    fitter.generateCovarianceMatrix(data_j_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::Correlated:
    fitter.generateCovarianceMatrix(data_j_inrange, CostType::Correlated);
    break;
  case CovarianceStrategy::CorrelatedBlock:
    fitter.generateCovarianceMatrix(data_j_ub_inrange, CostType::Correlated); //unbinned data
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
