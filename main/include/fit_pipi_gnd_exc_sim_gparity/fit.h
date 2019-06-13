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

  bool corr_mat_from_unbinned_data;
  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> const* corr_comb_j_unbinned;
  
  bool load_priors;
  std::string load_priors_file;

  MinimizerType minimizer;
  bool load_minimizer_params;
  std::string minimizer_params_file;

  bool load_bounds;
  std::string load_bounds_file;

  fitOptions(): load_frozen_fit_params(false),  write_covariance_matrix(false), 
		corr_mat_from_unbinned_data(false), corr_comb_j_unbinned(NULL), load_priors(false), 
		load_minimizer_params(false), load_bounds(false), minimizer(MinimizerType::MarquardtLevenberg){}
};

#define PRIOR_T_MEMBERS ( double, value )( double, weight )( int, param_idx )
struct PriorT{
  GENERATE_MEMBERS(PRIOR_T_MEMBERS);
  PriorT(): value(1.0), weight(0.2), param_idx(0){}
};
GENERATE_PARSER(PriorT, PRIOR_T_MEMBERS);

#define PRIOR_ARGS_MEMBERS ( std::vector<PriorT>, priors )
struct PriorArgs{
  GENERATE_MEMBERS(PRIOR_ARGS_MEMBERS);
  PriorArgs(): priors(1){}
};
GENERATE_PARSER(PriorArgs, PRIOR_ARGS_MEMBERS);

#define BOUND_ARGS_MEMBERS ( std::vector<boundedParameterTransform>, bounds )
struct BoundArgs{
  GENERATE_MEMBERS(BOUND_ARGS_MEMBERS);
  BoundArgs(): bounds(1){}
};
GENERATE_PARSER(BoundArgs, BOUND_ARGS_MEMBERS);



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


inline generalContainer getMinimizerParams(const fitOptions &opt){
  return getMinimizerParams(opt.minimizer, opt.load_minimizer_params, opt.minimizer_params_file);
}



//For testing bin size dependence it can be useful to fix the correlation matrix from the unbinned data and do a frozen fit
void generateFrozenCovMatFromUnbinnedData(simpleFitWrapper &fit,
					  const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j_unbinned,
					  const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j_binned){
  std::cout << "Generating frozen covariance matrix from binned jackknife sigma and unbinned jackknife correlation matrix" << std::endl;
  assert(corr_comb_j_unbinned.size() == corr_comb_j_binned.size());
  int ndata = corr_comb_j_binned.size();
  int nsample_binned = corr_comb_j_binned.value(0).size();
  NumericSquareMatrix<jackknifeDistributionD > corr(ndata);
  std::vector<jackknifeDistributionD > sigma(ndata);

  NumericSquareMatrix<double> cov_mn_unbinned(ndata);
  NumericSquareMatrix<double> cov_mn_binned(ndata);
  for(int i=0;i<ndata;i++)
    for(int j=i;j<ndata;j++){
      cov_mn_binned(i,j) = cov_mn_binned(j,i) = jackknifeDistribution<double>::covariance(corr_comb_j_binned.value(i), corr_comb_j_binned.value(j));
      cov_mn_unbinned(i,j) = cov_mn_unbinned(j,i) = jackknifeDistribution<double>::covariance(corr_comb_j_unbinned.value(i), corr_comb_j_unbinned.value(j));
    }

  NumericSquareMatrix<double> corr_mn_unbinned(ndata);
  for(int i=0;i<ndata;i++)
    for(int j=i;j<ndata;j++)
      corr_mn_unbinned(i,j) = corr_mn_unbinned(j,i) = cov_mn_unbinned(i,j)/sqrt( cov_mn_unbinned(i,i) * cov_mn_unbinned(j,j) );
  
  for(int i=0;i<ndata;i++){
    sigma[i] = jackknifeDistributionD(nsample_binned, sqrt(cov_mn_binned(i,i)));
    for(int j=i;j<ndata;j++)
      corr(i,j) = corr(j,i) = jackknifeDistributionD(nsample_binned, corr_mn_unbinned(i,j));
  }
  fit.importCorrelationMatrix(corr,sigma);
}
			
//In this version we allow the weights sigma to fluctuate between jackknife samples but the correlation matrix is fixed to that obtained from 
//the unbinned data. The idea is that this might capture much of the fluctuations that we would have with a true correlated fit
void generatePartiallyFrozenCovMatFromUnbinnedData(simpleFitWrapper &fit,
						   const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j_unbinned,
						   const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj_binned){
  std::cout << "Generating frozen covariance matrix from binned jackknife sigma and unbinned jackknife correlation matrix" << std::endl;
  assert(corr_comb_j_unbinned.size() == corr_comb_dj_binned.size());
  int ndata = corr_comb_dj_binned.size();
  int nsample_binned = corr_comb_dj_binned.value(0).size();
  NumericSquareMatrix<jackknifeDistributionD > corr(ndata);
  std::vector<jackknifeDistributionD > sigma(ndata);

  NumericSquareMatrix<double> cov_mn_unbinned(ndata);
  for(int i=0;i<ndata;i++)
    for(int j=i;j<ndata;j++)
      cov_mn_unbinned(i,j) = cov_mn_unbinned(j,i) = jackknifeDistributionD::covariance(corr_comb_j_unbinned.value(i), corr_comb_j_unbinned.value(j));
    
  NumericSquareMatrix<double> corr_mn_unbinned(ndata);
  for(int i=0;i<ndata;i++)
    for(int j=i;j<ndata;j++)
      corr_mn_unbinned(i,j) = corr_mn_unbinned(j,i) = cov_mn_unbinned(i,j)/sqrt( cov_mn_unbinned(i,i) * cov_mn_unbinned(j,j) );
  
  for(int i=0;i<ndata;i++){
    sigma[i] = sqrt( doubleJackknifeDistributionD::covariance(corr_comb_dj_binned.value(i),corr_comb_dj_binned.value(i) ) );
    for(int j=i;j<ndata;j++)
      corr(i,j) = corr(j,i) = jackknifeDistributionD(nsample_binned, corr_mn_unbinned(i,j));
  }
  fit.importCorrelationMatrix(corr,sigma);
}





//Note: nstate applies only for "MultiState" variants
void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 const correlationFunction<SimFitCoordGen,  blockDoubleJackknifeDistributionD> &corr_comb_bdj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int nstate, const int Lt, 
	 const int t_min, const int t_max,
	 const bool correlated, const CovarianceMatrix covariance_matrix,
	 const double Ascale, const double Cscale,
	 const fitOptions &opt = fitOptions()){
  
  std::unique_ptr<genericFitFuncBase> fitfunc = getFitFunc(ffunc, nstate, t_min, Lt, param_map.size(), Ascale, Cscale, params.sample(0));
  
  generalContainer min_params = getMinimizerParams(opt);

  simpleFitWrapper fit(*fitfunc, opt.minimizer, min_params);
  
  const int nsample = corr_comb_j.value(0).size();
 
  //Read frozen fit parameters
  if(opt.load_frozen_fit_params){
    std::cout << "Reading frozen fit params from " << opt.load_frozen_fit_params_file << std::endl;
    readFrozenParams(fit, opt.load_frozen_fit_params_file, nsample);
  }

  //Generate covariance matrix
  std::cout << "Generating and importing covariance matrix\n";
  CostType cost_type = correlated ? CostType::Correlated : CostType::Uncorrelated;  

  if(opt.corr_mat_from_unbinned_data){
    assert(covariance_matrix != CovarianceMatrix::BlockHybrid);
    //For testing what effect binning has on the covariance matrix separate from any underlying error dependence, here
    //we generate the covariance matrix using the correlation matrix from the unbinned data
    //If we use frozen fits the weights sigma will be computed from the binned jackknife data and fixed for all samples,
    //otherwise for unfrozen fits sigma will be computed from binned double-jackknife data. Note in both cases the same fixed correlation matrix is used
    assert(opt.corr_comb_j_unbinned != NULL);    

    if(covariance_matrix == CovarianceMatrix::Frozen) generateFrozenCovMatFromUnbinnedData(fit, *opt.corr_comb_j_unbinned, corr_comb_j);
    else generatePartiallyFrozenCovMatFromUnbinnedData(fit, *opt.corr_comb_j_unbinned, corr_comb_dj);
  }else{
    if(covariance_matrix == CovarianceMatrix::Frozen) fit.generateCovarianceMatrix(corr_comb_j, cost_type);
    else if(covariance_matrix == CovarianceMatrix::BlockHybrid) fit.generateCovarianceMatrix(corr_comb_dj, corr_comb_bdj, cost_type);
    else fit.generateCovarianceMatrix(corr_comb_dj, cost_type);
  }

  if(opt.write_covariance_matrix) fit.writeCovarianceMatrixHDF5(opt.write_covariance_matrix_file);
  
  //Add priors if necessary
  int Nprior;
  if(opt.load_priors){
    PriorArgs pargs;  parse(pargs, opt.load_priors_file);
    for(int p=0;p<pargs.priors.size();p++){
      std::cout << "Added prior value " << pargs.priors[p].value << " weight " 
		<<  pargs.priors[p].weight << " param " << pargs.priors[p].param_idx << std::endl;

	fit.addPrior(pargs.priors[p].value, pargs.priors[p].weight, pargs.priors[p].param_idx);
    }
    Nprior = pargs.priors.size();
  }

  //Add bounds if necessary
  if(opt.load_bounds){
    BoundArgs bargs;  parse(bargs, opt.load_bounds_file);
    for(int p=0;p<bargs.bounds.size();p++){
      std::cout << "Adding bound " << bargs.bounds[p].bound << " to param " << bargs.bounds[p].param 
		<< " with min " << bargs.bounds[p].min << " and/or max " <<bargs.bounds[p].max << std::endl; 

      fit.setBound(bargs.bounds[p]);
    }
  }

  //Run the fit
  std::cout << "Running fit routine\n";
  int dof;
  std::pair<jackknifeDistribution<double>, int> chisq_dof_nopriors;

  fit.fit(params, chisq, chisq_per_dof, dof, corr_comb_j, opt.load_priors ? &chisq_dof_nopriors : NULL);
  
  if(opt.load_priors) std::cout << "Chi^2 (excl priors) from jackknife = " << chisq_dof_nopriors.first << std::endl;    
}



void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 const correlationFunction<SimFitCoordGen,  blockDoubleJackknifeDistributionD> &corr_comb_bdj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int Lt, 
	 const int t_min, const int t_max,
	 const bool correlated, const CovarianceMatrix covariance_matrix,
	 const double Ascale, const double Cscale,
	 const fitOptions &opt = fitOptions()){
  assert(ffunc != FitFuncType::FSimGenMultiState);
  fit(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, corr_comb_bdj, ffunc, param_map, 0, Lt, t_min, t_max, correlated, covariance_matrix, Ascale, Cscale, opt);
}

CPSFIT_END_NAMESPACE

#endif
