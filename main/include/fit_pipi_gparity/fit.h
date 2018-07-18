#ifndef _FIT_PIPI_GPARITY_FIT_H_
#define _FIT_PIPI_GPARITY_FIT_H_

#include<fit.h>

typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunction;

struct pipiFitOptions{
  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  bool load_guess;
  std::string guess_file;

  pipiFitOptions(): load_frozen_fit_params(false), load_guess(false){}

  template<typename T>
  inline void import(const T &from){ 
#define I(A) A = from.A
    I(load_frozen_fit_params);
    I(load_frozen_fit_params_file);
    I(load_guess);
    I(guess_file);
#undef I
  }    
};


template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
jackknifeDistribution<typename FitFunc::Params> fit_corr_uncorr(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
								const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
								const int Lt, const int tsep_pipi, const double Ascale, const double Cscale,
								const pipiFitOptions &opt = pipiFitOptions()){
  typedef typename FitFunc::Params Params;
  Params guess;
  if(opt.load_guess){
    parse(guess, opt.guess_file);
    std::cout << "Loaded guess: " << guess << std::endl;
  }
  const int nsample = pipi_dj_vacsubbed_inrange.value(0).size();
  
  FitFunc fitfunc(Lt, tsep_pipi, Ascale, Cscale);

  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
  fitter<FitPolicies> fitter;
  fitter.importFitFunc(fitfunc);

  if(opt.load_frozen_fit_params)
    readFrozenParams(fitter, opt.load_frozen_fit_params_file, nsample);
  
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
								       const bool correlated, 
								       const int Lt, const int tsep_pipi, const double Ascale, const double Cscale,
								       const pipiFitOptions &opt = pipiFitOptions()){
  jackknifeDistribution<typename FitFunc::Params> params = correlated ? 
    fit_corr_uncorr<FitFunc,correlatedFitPolicy>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, Lt, tsep_pipi, Ascale, Cscale, opt) :
    fit_corr_uncorr<FitFunc,uncorrelatedFitPolicy>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, Lt, tsep_pipi, Ascale, Cscale, opt);

  jackknifeDistributionD Epipi(params.size(), [&](const int s){ return params.sample(s).pipiEnergy(); });
  jackknifeDistributionD constant(params.size(), [&](const int s){ return params.sample(s).constant(); });
  return std::pair<jackknifeDistributionD,jackknifeDistributionD>(std::move(Epipi),std::move(constant));
}

//returns a pair containing the pipi energy and the constant term
inline std::pair<jackknifeDistributionD,jackknifeDistributionD> fit(const jackknifeCorrelationFunction &pipi_j_vacsubbed_inrange,
								    const doubleJackCorrelationFunction &pipi_dj_vacsubbed_inrange,
								    const PiPiFitFunction fitfunc, const bool correlated,
								    const int Lt, const int tsep_pipi, const double Ascale, const double Cscale,
								    const pipiFitOptions &opt = pipiFitOptions()){
  switch(fitfunc){
  case PiPiFitFunction::FCoshPlusConstant:
    return fit_ff<FitCoshPlusConstant>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, correlated, Lt, tsep_pipi, Ascale, Cscale, opt);
  case PiPiFitFunction::FCoshPlusConstantDoubleExp:
    return fit_ff<FitCoshPlusConstantDoubleExp>(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange, correlated, Lt, tsep_pipi, Ascale, Cscale, opt);
  default:
    error_exit(std::cout << "Unknown fitfunc " << fitfunc << std::endl);
  }
};

#endif
