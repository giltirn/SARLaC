#ifndef _PIPI_GND_EXC_SIM_FIT_FIT_BOOTSTRAP_H
#define _PIPI_GND_EXC_SIM_FIT_FIT_BOOTSTRAP_H

#include<config.h>
#include<utils/macros.h>
#include<fit_pipi_gnd_exc_sigma_sim_gparity/fit.h>

CPSFIT_START_NAMESPACE

//Note: nstate applies only for "MultiState" variants
void fit(bootstrapDistribution<taggedValueContainer<double,std::string> > &params, 
	 bootstrapDistributionD &chisq, bootstrapDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  bootstrapDistributionD> &corr_comb_b,
	 const correlationFunction<SimFitCoordGen,  bootJackknifeDistributionD> &corr_comb_bj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int nstate, const int Lt, 
	 const int t_min, const int t_max,
	 const bool correlated, const CovarianceMatrix covariance_matrix,
	 const double Ascale, const double Cscale,
	 const fitOptions &opt = fitOptions()){
  
  assert(!opt.corr_mat_from_unbinned_data);

  std::unique_ptr<genericFitFuncBase> fitfunc = getFitFunc(ffunc, nstate, t_min, Lt, param_map.size(), Ascale, Cscale, params.best());
  
  generalContainer min_params = getMinimizerParams(opt);

  simpleFitWrapper<bootstrapDistributionD> fit(*fitfunc, opt.minimizer, min_params);
  
  const int nboot = corr_comb_b.value(0).size();
 
  //Read frozen fit parameters
  if(opt.load_frozen_fit_params){
    std::cout << "Reading frozen fit params from " << opt.load_frozen_fit_params_file << std::endl;
    readFrozenParams(fit, opt.load_frozen_fit_params_file, nboot);
  }

  //Generate covariance matrix
  std::cout << "Generating and importing covariance matrix\n";
  CostType cost_type = correlated ? CostType::Correlated : CostType::Uncorrelated;  

  if(covariance_matrix == CovarianceMatrix::Frozen) fit.generateCovarianceMatrix(corr_comb_b, cost_type);
  else if(covariance_matrix == CovarianceMatrix::Regular) fit.generateCovarianceMatrix(corr_comb_bj, cost_type);
  else assert(0);

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
  std::pair<bootstrapDistributionD, int> chisq_dof_nopriors(bootstrapDistributionD(bootstrapInitType(nboot)), 0);

  fit.fit(params, chisq, chisq_per_dof, dof, corr_comb_b, opt.load_priors ? &chisq_dof_nopriors : NULL);
  
  if(opt.load_priors) std::cout << "Chi^2 (excl priors) from bootstrap = " << chisq_dof_nopriors.first << std::endl;    
}



CPSFIT_END_NAMESPACE

#endif
