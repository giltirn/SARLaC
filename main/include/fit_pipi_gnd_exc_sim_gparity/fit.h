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

  fitOptions(): load_frozen_fit_params(false),  write_covariance_matrix(false), load_priors(false){}
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

    if(opt.load_priors){
      PriorArgs pargs;  parse(pargs, opt.load_priors_file);
      for(int p=0;p<pargs.priors.size();p++){
	std::cout << "Added prior value " << pargs.priors[p].value << " weight " <<  pargs.priors[p].weight << " param " << pargs.priors[p].param_idx << std::endl;
	fit.addPrior(pargs.priors[p].value, pargs.priors[p].weight, pargs.priors[p].param_idx);
      }
    }

    std::cout << "Running fit routine\n";
    fit.fit(params, chisq, chisq_per_dof, corr_comb_j);
}  

void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int Lt, const double Ascale, const double Cscale,
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
  }else{
    assert(0);
  }
}


CPSFIT_END_NAMESPACE

#endif
