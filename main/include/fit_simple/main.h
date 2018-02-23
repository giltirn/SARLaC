#ifndef _FIT_SIMPLE_MAIN_H_
#define _FIT_SIMPLE_MAIN_H_

template<typename FitFunc, template<typename> class CostFunctionPolicy, typename ArgsType, typename CMDlineType>
void fitSpecFFcorr(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline){
  typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, CostFunctionPolicy>::type FitPolicies;
  typedef fitter<FitPolicies> Fitter;
  typedef typename FitPolicies::FitParameterDistribution FitParameterDistribution;
  
  const int nt_fit = args.t_max - args.t_min + 1;

  doubleJackknifeCorrelationFunctionD data_dj_inrange(nt_fit, 
						      [&](const int i){ return typename doubleJackknifeCorrelationFunctionD::ElementType(args.t_min + i, data_dj.value(args.t_min + i)); }
						      );
  
  jackknifeCorrelationFunctionD data_j_inrange(nt_fit, 
					       [&](const int i){ return typename jackknifeCorrelationFunctionD::ElementType(args.t_min + i, data_j.value(args.t_min + i)); }
					       );

  FitFunc* fitfunc = FitFuncPolicy<FitFunc,ArgsType>::get(args);
  Fitter fit;
  fit.importFitFunc(*fitfunc);

  importCostFunctionParameters<CostFunctionPolicy,FitPolicies> prepare(fit,data_dj_inrange);

  typename FitFunc::ParameterType guess = fitfunc->guess();
  if(cmdline.load_guess){
    parse(guess,cmdline.guess_file);
  }

  const int nsample = data_j.value(0).size();
  FitParameterDistribution params(nsample, guess);
  jackknifeDistributionD chisq(nsample);
  jackknifeDistributionD chisq_per_dof(nsample);

  fit.fit(params, chisq, chisq_per_dof, data_j_inrange);

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
#ifdef HAVE_HDF5
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(params, "params.hdf5");
#endif

  FitFuncPolicy<FitFunc,ArgsType>::plot(args,*fitfunc,data_j,params);
  delete fitfunc;
}

#endif
