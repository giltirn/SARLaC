#ifndef _FIT_PIPI_GPARITY_FIT_H_
#define _FIT_PIPI_GPARITY_FIT_H_

#include<fit.h>

typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunction;

template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
jackknifeDistribution<typename FitFunc::Params> fit_corr_uncorr(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
								const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
								const Args &args, const CMDline &cmdline){

  typedef typename FitFunc::Params Params;
  Params guess;
  if(cmdline.load_guess){
    parse(guess, cmdline.guess_file);
    std::cout << "Loaded guess: " << guess << std::endl;
  }
  const int nsample = pipi_dj_vacsubbed_inrange.value(0).size();
  
  FitFunc fitfunc(args.Lt, args.tsep_pipi, args.Ascale, args.Cscale);

  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
  fitter<FitPolicies> fitter;
  fitter.importFitFunc(fitfunc);

  if(cmdline.load_frozen_fit_params)
    readFrozenParams(fitter, cmdline.load_frozen_fit_params_file, nsample);
  
  importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, pipi_dj_vacsubbed_inrange);
    
  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq;
  jackknifeDistributionD chisq_per_dof;
  fitter.fit(params, chisq, chisq_per_dof, pipi_j_vacsubbed_inrange);

  double dof = chisq.sample(0)/chisq_per_dof.sample(0);
  
  jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });
  
  distributionPrint<jackknifeDistribution<Params> >::printer(new pipiParamsPrinter<FitFunc>);

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "Dof: " << dof << std::endl;
  std::cout << "P-value: " << pvalue << std::endl;

#ifdef HAVE_HDF5
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(pvalue, "pvalue.hdf5");
#endif
  
  return params;
}
  
template<typename FitFunc>
inline std::pair<jackknifeDistributionD,jackknifeDistributionD> fit_ff(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
								       const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
								       const Args &args, const CMDline &cmdline){
  jackknifeDistribution<typename FitFunc::Params> params = args.correlated ? 
    fit_corr_uncorr<FitFunc,correlatedFitPolicy>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, args, cmdline) :
    fit_corr_uncorr<FitFunc,uncorrelatedFitPolicy>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, args, cmdline);

  jackknifeDistributionD Epipi(params.size(), [&](const int s){ return params.sample(s).pipiEnergy(); });
  jackknifeDistributionD constant(params.size(), [&](const int s){ return params.sample(s).constant(); });
  return std::pair<jackknifeDistributionD,jackknifeDistributionD>(std::move(Epipi),std::move(constant));
}

//returns a pair containing the pipi energy and the constant term
inline std::pair<jackknifeDistributionD,jackknifeDistributionD> fit(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
		const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
		const Args &args, const CMDline &cmdline){
  switch(args.fitfunc){
  case FCoshPlusConstant:
    return fit_ff<FitCoshPlusConstant>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange,args,cmdline);
  case FCoshPlusConstantDoubleExp:
    return fit_ff<FitCoshPlusConstantDoubleExp>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange,args,cmdline);
  default:
    error_exit(std::cout << "Unknown fitfunc " << args.fitfunc << std::endl);
  }
};

#endif
