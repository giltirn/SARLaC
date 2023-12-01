#ifndef _PIPI_GND_EXC_SIM_FIT_FIT_CENTRAL_H
#define _PIPI_GND_EXC_SIM_FIT_FIT_CENTRAL_H

#include<config.h>
#include<utils/macros.h>

#include "fit.h"

SARLAC_START_NAMESPACE


struct fitCentralOptions{
  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;
 
  bool load_priors;
  std::string load_priors_file;

  MinimizerType minimizer;
  bool load_minimizer_params;
  std::string minimizer_params_file;

  bool load_bounds;
  std::string load_bounds_file;

  fitCentralOptions(): load_frozen_fit_params(false),
		       load_priors(false), load_minimizer_params(false), load_bounds(false), 
		       minimizer(MinimizerType::MarquardtLevenberg){}
};

//Fit just the central value. Jackknife data is used to obtain the covariance matrix
//Note as we do not compute an error, binning is irrelevant
void fitCentral(taggedValueContainer<double,std::string> &params, double &chisq, double &chisq_per_dof, 
		const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
		const std::unique_ptr<genericFitFuncBase> &fitfunc, const std::unordered_map<std::string,size_t> &param_map,
		const int nstate, const int Lt, 
		const int t_min, const int t_max,
		const bool correlated,
		const double Ascale, const double Cscale,
		const fitCentralOptions &opt = fitCentralOptions()){
  
  generalContainer min_params = getMinimizerParams(opt.minimizer, opt.load_minimizer_params, opt.minimizer_params_file);

  simpleSingleFitWrapper fit(*fitfunc, opt.minimizer, min_params);
  
  const int nsample = corr_comb_j.value(0).size();
 
  //Read frozen fit parameters
  if(opt.load_frozen_fit_params){
    std::cout << "Reading frozen fit params from " << opt.load_frozen_fit_params_file << std::endl;
    readFrozenParams(fit, opt.load_frozen_fit_params_file, nsample);
  }

  //Generate covariance matrix (only correlated or uncorrelated are options)
  std::cout << "Generating and importing covariance matrix\n";
  CostType cost_type = correlated ? CostType::Correlated : CostType::Uncorrelated;  
  
  fit.generateCovarianceMatrix(corr_comb_j, cost_type);
  
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

  
  correlationFunction<SimFitCoordGen,  double> corr_comb_cen(corr_comb_j.size());
  for(int i=0;i<corr_comb_j.size();i++){
    corr_comb_cen.coord(i) = corr_comb_j.coord(i);
    corr_comb_cen.value(i) = corr_comb_j.value(i).mean();
  }
  
  fit.fit(params, chisq, chisq_per_dof, dof, corr_comb_cen);
}

void fitCentral(taggedValueContainer<double,std::string> &params, double &chisq, double &chisq_per_dof, 
		const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
		FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
		const int nstate, const int Lt, 
		const int t_min, const int t_max,
		const bool correlated,
		const double Ascale, const double Cscale,
		const fitCentralOptions &opt = fitCentralOptions()){
  
  std::unique_ptr<genericFitFuncBase> fitfunc = getFitFunc(ffunc, nstate, t_min, Lt, param_map.size(), Ascale, Cscale, params);
  
  fitCentral(params, chisq, chisq_per_dof, corr_comb_j, fitfunc, param_map, nstate, Lt, t_min, t_max, correlated, Ascale, Cscale, opt);
}


SARLAC_END_NAMESPACE

#endif
