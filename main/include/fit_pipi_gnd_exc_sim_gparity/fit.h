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

  MinimizerType minimizer;
  bool load_minimizer_params;
  std::string minimizer_params_file;

  fitOptions(): load_frozen_fit_params(false),  write_covariance_matrix(false), load_priors(false), 
		load_minimizer_params(false), minimizer(MinimizerType::MarquardtLevenberg){}
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

std::unique_ptr<genericFitFuncBase> getFitFunc(const FitFuncType type, const int nstate, const int t_min, const int Lt, 
					       const int nparam, const double Ascale, const double Cscale,
					       const taggedValueContainer<double,std::string> &psetup){
  typedef std::unique_ptr<genericFitFuncBase> PtrType;

#define FDEF_BASIC(ENUM, CLASS) \
  case(FitFuncType:: ENUM):						\
    return PtrType(new genericFitFuncWrapper<CLASS>(CLASS(Lt, nparam, Ascale, Cscale), psetup)); break
  
#define FDEF_NSTATE(ENUM, CLASS)					\
  case(FitFuncType:: ENUM):						\
    return PtrType(new genericFitFuncWrapper<CLASS>(CLASS(nstate, Lt, nparam, Ascale, Cscale), psetup)); break
  
  switch(type){
    FDEF_BASIC(FSimGenOneState, FitSimGenOneState);
    FDEF_BASIC(FSimGenTwoState, FitSimGenTwoState);
    FDEF_BASIC(FSimGenThreeState, FitSimGenThreeState);
    FDEF_BASIC(FSimGenThreeStateLogEdiff, FitSimGenThreeStateLogEdiff);
    FDEF_NSTATE(FSimGenMultiState, FitSimGenMultiState);
    FDEF_NSTATE(FSimGenMultiStateLogEdiff, FitSimGenMultiStateLogEdiff);
    FDEF_NSTATE(FSimGenMultiStateCparam, FitSimGenMultiStateCparam);
    FDEF_NSTATE(FSimGenMultiStateSub, FitSimGenMultiStateSub);    
    case(FitFuncType::FSimGenMultiStateTminSub):
    case(FitFuncType::FSimGenMultiStateTminSubForceZero):
      return PtrType(new genericFitFuncWrapper<FitSimGenMultiStateTminSub>(FitSimGenMultiStateTminSub(nstate, Lt, nparam, t_min, Ascale, Cscale), psetup)); break;  
    default:
      assert(0);
  }

}

template<typename T>
generalContainer getMinimizerParamsT(const fitOptions &opt){
  generalContainer min_params;
  T mp; mp.verbose = true;
  if(opt.load_minimizer_params){
    parse(mp, opt.minimizer_params_file);
    std::cout << "Loaded minimizer params: " << mp << std::endl;
  }
  min_params = mp;
  return min_params;
}
generalContainer getMinimizerParams(const fitOptions &opt){
  switch(opt.minimizer){
  case MinimizerType::MarquardtLevenberg:
    return getMinimizerParamsT<MarquardtLevenbergParameters<double> >(opt);
  case MinimizerType::GSLtrs:
    return getMinimizerParamsT<GSLtrsMinimizerParams>(opt);
  case MinimizerType::GSLmultimin:
    return getMinimizerParamsT<GSLmultidimMinimizerParams>(opt);
  case MinimizerType::Minuit2:
#ifdef HAVE_MINUIT2
    return getMinimizerParamsT<Minuit2minimizerParams>(opt);
#else
    error_exit(std::cout << "Library not compiled with Minuit2\n");
#endif
  default:
    assert(0);
  }
}


//Note: nstate applies only for "MultiState" variants
void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int nstate, const int Lt, 
	 const int t_min, const int t_max,
	 const bool correlated,
	 const double Ascale, const double Cscale,
	 const fitOptions &opt = fitOptions()){
  
  std::unique_ptr<genericFitFuncBase> fitfunc = getFitFunc(ffunc, nstate, t_min, Lt, param_map.size(), Ascale, Cscale, params.sample(0));
  
  generalContainer min_params = getMinimizerParams(opt);

  simpleFitWrapper fit(*fitfunc, opt.minimizer, min_params);
  
  const int nsample = corr_comb_j.value(0).size();
 
  if(opt.load_frozen_fit_params){
    std::cout << "Reading frozen fit params from " << opt.load_frozen_fit_params_file << std::endl;
    readFrozenParams(fit, opt.load_frozen_fit_params_file, nsample);
  }

  std::cout << "Generating and importing covariance matrix\n";
  CostType cost_type = correlated ? CostType::Correlated : CostType::Uncorrelated;  
  fit.generateCovarianceMatrix(corr_comb_dj, cost_type);

  if(opt.write_covariance_matrix) fit.writeCovarianceMatrixHDF5(opt.write_covariance_matrix_file);

  int Nprior;
  if(opt.load_priors){
    PriorArgs pargs;  parse(pargs, opt.load_priors_file);
    for(int p=0;p<pargs.priors.size();p++){
      std::cout << "Added prior value " << pargs.priors[p].value << " weight " <<  pargs.priors[p].weight << " param " << pargs.priors[p].param_idx << std::endl;
	fit.addPrior(pargs.priors[p].value, pargs.priors[p].weight, pargs.priors[p].param_idx);
    }
    Nprior = pargs.priors.size();
  }

  std::cout << "Running fit routine\n";
  int dof;
  std::pair<jackknifeDistribution<double>, int> chisq_dof_nopriors;

  fit.fit(params, chisq, chisq_per_dof, dof, corr_comb_j, opt.load_priors ? &chisq_dof_nopriors : NULL);
  
  if(opt.load_priors) std::cout << "Chi^2 (excl priors) from jackknife = " << chisq_dof_nopriors.first << std::endl;    
}



void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int Lt, 
	 const int t_min, const int t_max,
	 const bool correlated, 
	 const double Ascale, const double Cscale,
	 const fitOptions &opt = fitOptions()){
  assert(ffunc != FitFuncType::FSimGenMultiState);
  fit(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, ffunc, param_map, 0, Lt, t_min, t_max, correlated, Ascale, Cscale, opt);
}

CPSFIT_END_NAMESPACE

#endif
