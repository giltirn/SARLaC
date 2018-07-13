#ifndef _SIGMA_FIT_H___
#define _SIGMA_FIT_H___

GENERATE_ENUM_AND_PARSER(SigmaFitFunction, (FCoshPlusConstant)(FCoshPlusConstantDoubleExp) );

struct SigmaFitArgs{
  bool load_guess;
  std::string guess_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  double Ascale;
  double Cscale;

  int Lt;

  bool correlated;
  
  SigmaFitFunction fitfunc;

  SigmaFitArgs(): load_guess(false), load_frozen_fit_params(false), Ascale(1e13), Cscale(1e13), Lt(64), correlated(true), fitfunc(SigmaFitFunction::FCoshPlusConstant){}
};

template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
jackknifeDistribution<typename FitFunc::Params> fit_corr_uncorr_sigma(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
								      const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
								      const SigmaFitArgs &args){

  typedef typename FitFunc::Params Params;
  Params guess;
  if(args.load_guess){
    parse(guess, args.guess_file);
    std::cout << "Loaded guess: " << guess << std::endl;
  }
  const int nsample = pipi_dj_vacsubbed_inrange.value(0).size();
  
  FitFunc fitfunc(args.Lt, 0, args.Ascale, args.Cscale); //use same fitfuncs as pipi, just set tsep_pipi to 0

  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
  fitter<FitPolicies> fitter;
  fitter.importFitFunc(fitfunc);

  if(args.load_frozen_fit_params)
    readFrozenParams(fitter, args.load_frozen_fit_params_file, nsample);
  
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
inline std::pair<jackknifeDistributionD,jackknifeDistributionD> fit_ff_sigma(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
								       const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
								       const  SigmaFitArgs &args){
  jackknifeDistribution<typename FitFunc::Params> params = args.correlated ? 
    fit_corr_uncorr_sigma<FitFunc,correlatedFitPolicy>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, args) :
    fit_corr_uncorr_sigma<FitFunc,uncorrelatedFitPolicy>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, args);

  jackknifeDistributionD Epipi(params.size(), [&](const int s){ return params.sample(s).pipiEnergy(); });
  jackknifeDistributionD constant(params.size(), [&](const int s){ return params.sample(s).constant(); });
  return std::pair<jackknifeDistributionD,jackknifeDistributionD>(std::move(Epipi),std::move(constant));
}

//returns a pair containing the pipi energy and the constant term
inline std::pair<jackknifeDistributionD,jackknifeDistributionD> fit_sigma(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
									  const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
									  const SigmaFitArgs &args){
  switch(args.fitfunc){
  case SigmaFitFunction::FCoshPlusConstant:
    return fit_ff_sigma<FitCoshPlusConstant>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange,args);
  case SigmaFitFunction::FCoshPlusConstantDoubleExp:
    return fit_ff_sigma<FitCoshPlusConstantDoubleExp>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange,args);
  default:
    error_exit(std::cout << "Unknown fitfunc " << args.fitfunc << std::endl);
  }
};

#endif
