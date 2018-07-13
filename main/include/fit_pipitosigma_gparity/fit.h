#ifndef _PIPI_TO_SIGMA_FIT_H_
#define _PIPI_TO_SIGMA_FIT_H_

GENERATE_ENUM_AND_PARSER(PiPiToSigmaFitFunction, (FCoshPlusConstant) );

struct PiPiToSigmaFitArgs{
  bool load_guess;
  std::string guess_file;

  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  double Ascale;
  double Cscale;

  int Lt;

  int tsep_pipi;

  bool correlated;
  
  PiPiToSigmaFitFunction fitfunc;

  PiPiToSigmaFitArgs(): load_guess(false), load_frozen_fit_params(false), Ascale(1e13), Cscale(1e13), Lt(64), correlated(true), fitfunc(PiPiToSigmaFitFunction::FCoshPlusConstant), tsep_pipi(4){}
};


template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
jackknifeDistribution<typename FitFunc::Params> fit_corr_uncorr_pipi_to_sigma(const jackknifeCorrelationFunction &data_j_vacsubbed_inrange,
								      const doubleJackCorrelationFunction &data_dj_vacsubbed_inrange,
								      const PiPiToSigmaFitArgs &args){

  typedef typename FitFunc::Params Params;
  Params guess;
  if(args.load_guess){
    parse(guess, args.guess_file);
    std::cout << "Loaded guess: " << guess << std::endl;
  }
  const int nsample = data_dj_vacsubbed_inrange.value(0).size();
  
  FitFunc fitfunc(args.Lt, args.tsep_pipi/2, args.Ascale, args.Cscale); //use same fitfuncs as pipi, just set tsep_pipi to half its value to account for there being only one pipi operator

  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
  fitter<FitPolicies> fitter;
  fitter.importFitFunc(fitfunc);

  if(args.load_frozen_fit_params)
    readFrozenParams(fitter, args.load_frozen_fit_params_file, nsample);
  
  importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, data_dj_vacsubbed_inrange);
    
  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq;
  jackknifeDistributionD chisq_per_dof;
  fitter.fit(params, chisq, chisq_per_dof, data_j_vacsubbed_inrange);

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
inline std::pair<jackknifeDistributionD,jackknifeDistributionD> fit_ff_pipi_to_sigma(const jackknifeCorrelationFunction &data_j_vacsubbed_inrange,
								       const doubleJackCorrelationFunction &data_dj_vacsubbed_inrange,
								       const PiPiToSigmaFitArgs &args){
  jackknifeDistribution<typename FitFunc::Params> params = args.correlated ? 
    fit_corr_uncorr_pipi_to_sigma<FitFunc,correlatedFitPolicy>(data_j_vacsubbed_inrange,data_dj_vacsubbed_inrange, args) :
    fit_corr_uncorr_pipi_to_sigma<FitFunc,uncorrelatedFitPolicy>(data_j_vacsubbed_inrange,data_dj_vacsubbed_inrange, args);

  jackknifeDistributionD Epipi(params.size(), [&](const int s){ return params.sample(s).pipiEnergy(); });
  jackknifeDistributionD constant(params.size(), [&](const int s){ return params.sample(s).constant(); });
  return std::pair<jackknifeDistributionD,jackknifeDistributionD>(std::move(Epipi),std::move(constant));
}

//returns a pair containing the pipi energy and the constant term
inline std::pair<jackknifeDistributionD,jackknifeDistributionD> fit_pipi_to_sigma(const jackknifeCorrelationFunction &data_j_vacsubbed_inrange,
										  const doubleJackCorrelationFunction &data_dj_vacsubbed_inrange,
										  const PiPiToSigmaFitArgs &args){
  switch(args.fitfunc){
  case PiPiToSigmaFitFunction::FCoshPlusConstant:
    return fit_ff_pipi_to_sigma<FitCoshPlusConstant>(data_j_vacsubbed_inrange,data_dj_vacsubbed_inrange,args);
  default:
    error_exit(std::cout << "Unknown fitfunc " << args.fitfunc << std::endl);
  }
};



#endif
