#ifndef _PIPI_GND_EXC_SIM_FIT_FIT_H
#define _PIPI_GND_EXC_SIM_FIT_FIT_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

struct fitOptions{
  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  bool write_covariance_matrix;
  std::string write_covariance_matrix_file;

  bool load_priors;
  std::string load_priors_file;

  bool load_mlparams;
  std::string mlparams_file;

  fitOptions(): load_frozen_fit_params(false),  write_covariance_matrix(false), load_priors(false), load_mlparams(false){}
};

#define PRIOR_T_MEMBERS ( double, value )( double, weight )( int, param_idx )
struct PriorT{
  GENERATE_MEMBERS(PRIOR_T_MEMBERS)
  PriorT(): value(1.0), weight(0.2), param_idx(0){}
};
GENERATE_PARSER(PriorT, PRIOR_T_MEMBERS)

#define PRIOR_ARGS_MEMBERS ( std::vector<PriorT>, priors )
struct PriorArgs{
  GENERATE_MEMBERS(PRIOR_ARGS_MEMBERS)
  PriorArgs(): priors(1){}
};
GENERATE_PARSER(PriorArgs, PRIOR_ARGS_MEMBERS)



template<typename FitFunc>
void fit_ff(jackknifeDistribution<typename FitFunc::Params> &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	    const correlationFunction<SimFitCoordGen, jackknifeDistributionD> &corr_comb_j,
	    const correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> &corr_comb_dj,
	    const FitFunc &fitfunc, const fitOptions &opt = fitOptions()){
    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
    const int nsample = corr_comb_j.value(0).size();

    fitter<FitPolicies> fit;
    fit.importFitFunc(fitfunc);
    
    if(opt.load_frozen_fit_params){
      std::cout << "Reading frozen fit params from " << opt.load_frozen_fit_params_file << std::endl;
      readFrozenParams(fit, opt.load_frozen_fit_params_file, nsample, &params.sample(0));
    }else{
      std::cout << "Inserting default param instance to frozen fitter" << std::endl;
      fit.freeze({},params); //parameter type size is dynamic, thus for the internals to properly construct a full parameter vector it needs something to copy from even if no freezing is performed
    }

    std::cout << "Generating and importing covariance matrix\n";
    importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import(fit, corr_comb_dj);

    if(opt.write_covariance_matrix) import.writeCovarianceMatrixHDF5(opt.write_covariance_matrix_file);

    int Nprior;
    if(opt.load_priors){
      PriorArgs pargs;  parse(pargs, opt.load_priors_file);
      for(int p=0;p<pargs.priors.size();p++){
	std::cout << "Added prior value " << pargs.priors[p].value << " weight " <<  pargs.priors[p].weight << " param " << pargs.priors[p].param_idx << std::endl;
	fit.addPrior(pargs.priors[p].value, pargs.priors[p].weight, pargs.priors[p].param_idx);
      }
      Nprior = pargs.priors.size();
    }
    if(opt.load_mlparams){
      typename fitter<FitPolicies>::minimizerParamsType minparams;
      parse(minparams, opt.mlparams_file);
      std::cout << "Loaded minimizer params: " << minparams << std::endl;
      fit.setMinimizerParams(minparams);
    }


    std::cout << "Running fit routine\n";
    fit.fit(params, chisq, chisq_per_dof, corr_comb_j);

    if(opt.load_priors){
      //Compute chi^2 and chi^2/dof without priors
      typedef correlationFunction<SimFitCoordGen, double> cenCorrFunc;
      int N = corr_comb_j.size();

      {
	cenCorrFunc data_cen(N, [&](const int i){ return cenCorrFunc::ElementType(corr_comb_j.coord(i), corr_comb_j.value(i).best()); });
      
	NumericSquareMatrix<double> inv_corr(N);
	std::vector<double> sigma(N);
	for(int i=0;i<N;i++){
	  sigma[i] = import.sigma[i].best();

	  for(int j=0;j<N;j++)
	    inv_corr(i,j) = import.inv_corr(i,j).best();
	}
	CorrelatedChisqCostFunction<FitFunc, cenCorrFunc> costfunc(fitfunc, data_cen, sigma, inv_corr);
	double chisq_noprior = costfunc.cost(params.best());
	double chisq_per_dof_noprior = chisq_noprior/(costfunc.Ndof() + Nprior);

	std::cout << "Chi^2 (excl priors) = " << chisq_noprior << std::endl;
	std::cout << "Chi^2/dof (excl priors) = " << chisq_per_dof_noprior << std::endl;
      }


      int nsample = params.size();
      jackknifeDistributionD chisq_noprior_j(nsample);
      for(int s=0;s<nsample;s++){	
	cenCorrFunc data_j(N, [&](const int i){ return cenCorrFunc::ElementType(corr_comb_j.coord(i), corr_comb_j.value(i).sample(s)); });
      
	NumericSquareMatrix<double> inv_corr(N);
	std::vector<double> sigma(N);
	for(int i=0;i<N;i++){
	  sigma[i] = import.sigma[i].sample(s);

	  for(int j=0;j<N;j++)
	    inv_corr(i,j) = import.inv_corr(i,j).sample(s);
	}
	CorrelatedChisqCostFunction<FitFunc, cenCorrFunc> costfunc(fitfunc, data_j, sigma, inv_corr);
	chisq_noprior_j.sample(s) = costfunc.cost(params.sample(s));
      }
      std::cout << "Chi^2 (excl priors) from jackknife = " << chisq_noprior_j << std::endl;

    }
}

//Note: nstate applies only for "MultiState" variants
void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int nstate, const int Lt, const double Ascale, const double Cscale,
	 const fitOptions &opt = fitOptions()){

  if(ffunc == FitFuncType::FSimGenOneState){
    typedef FitSimGenOneState FitFunc;
    FitFunc fitfunc(Lt, param_map.size(), Ascale, Cscale);
    return fit_ff<FitFunc>(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, fitfunc, opt);
  }else if(ffunc == FitFuncType::FSimGenTwoState){
    typedef FitSimGenTwoState FitFunc;
    FitFunc fitfunc(Lt, param_map.size(), Ascale, Cscale);
    return fit_ff<FitFunc>(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, fitfunc, opt);
  }else if(ffunc == FitFuncType::FSimGenThreeState){
    typedef FitSimGenThreeState FitFunc;
    FitFunc fitfunc(Lt, param_map.size(), Ascale, Cscale);
    return fit_ff<FitFunc>(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, fitfunc, opt);
  }else if(ffunc == FitFuncType::FSimGenThreeStateLogEdiff){
    typedef FitSimGenThreeStateLogEdiff FitFunc;
    FitFunc fitfunc(Lt, param_map.size(), Ascale, Cscale);
    return fit_ff<FitFunc>(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, fitfunc, opt);
  }else if(ffunc == FitFuncType::FSimGenMultiState){
    typedef FitSimGenMultiState FitFunc;
    FitFunc fitfunc(nstate, Lt, param_map.size(), Ascale, Cscale);
    return fit_ff<FitFunc>(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, fitfunc, opt);
  }else{
    assert(0);
  }
}
void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int Lt, const double Ascale, const double Cscale,
	 const fitOptions &opt = fitOptions()){
  assert(ffunc != FitFuncType::FSimGenMultiState);
  fit(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, ffunc, param_map, 0, Lt, Ascale, Cscale, opt);
}

CPSFIT_END_NAMESPACE

#endif
