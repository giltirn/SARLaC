#ifndef _FIT_SIMPLE_MAIN_H_
#define _FIT_SIMPLE_MAIN_H_

template<typename DistributionType>
correlationFunction<double, DistributionType> resampleAndCombine(const std::vector<rawDataCorrelationFunctionD> &channels_raw, const Args &args, const CMDline &cmdline){
  const int nchannel = channels_raw.size();

  std::vector<correlationFunction<double, DistributionType> > channels_r(nchannel);
  for(int i=0;i<nchannel;i++){
    channels_r[i] = correlationFunction<double, DistributionType>(args.Lt, [&](const int t){
	return typename correlationFunction<double, DistributionType>::ElementType(t, DistributionType(channels_raw[i].value(t).bin(args.bin_size)));								  });
  }
  correlationFunction<double, DistributionType> out(args.Lt);

  applyCombination(out,channels_r,args.combination);
  applyTimeDep(out, args.outer_time_dep, args.Lt);

  return out;
}

template<typename DataSeriesType>
inline DataSeriesType getDataInRange(const DataSeriesType &data, const int tmin, const int tmax){
  return DataSeriesType(tmax - tmin + 1, [&](const int i){ return data[tmin + i]; });
}
  
template<typename FitFunc, template<typename> class CostFunctionPolicy, typename ArgsType, typename CMDlineType>
void fitSpecFFcorr(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const ArgsType &args, const CMDlineType &cmdline){
  typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, CostFunctionPolicy>::type FitPolicies;
  typedef fitter<FitPolicies> Fitter;
  typedef typename FitPolicies::FitParameterDistribution FitParameterDistribution;
  
  jackknifeCorrelationFunctionD data_j_inrange = getDataInRange(data_j, args.t_min, args.t_max);
  doubleJackknifeCorrelationFunctionD data_dj_inrange;
  if(args.covariance_strategy != CovarianceStrategy::FrozenCorrelated) data_dj_inrange = getDataInRange(data_dj, args.t_min, args.t_max);

  std::cout << "All data\n";
  for(int i=0;i<data_j.size();i++){
    std::cout << (int)data_j.coord(i) << " " << data_j.value(i) << std::endl;
  }

  std::cout << "Data in range\n";
  for(int i=0;i<data_j_inrange.size();i++){
    std::cout << (int)data_j_inrange.coord(i) << " " << data_j_inrange.value(i) << std::endl;
  }

  FitFunc* fitfunc = FitFuncPolicy<FitFunc,ArgsType>::get(args);
  Fitter fit;

  if(cmdline.load_mlparams){
    typename Fitter::minimizerParamsType minparams;
    parse(minparams, cmdline.mlparams_file);
    std::cout << "Loaded minimizer params: " << minparams << std::endl;
    fit.setMinimizerParams(minparams);
  }

  fit.importFitFunc(*fitfunc);

  importCostFunctionParameters<CostFunctionPolicy,FitPolicies> prepare(fit,data_j_inrange, data_dj_inrange);

  typename FitFunc::ParameterType guess = fitfunc->guess();
  if(cmdline.load_guess){
    parse(guess,cmdline.guess_file);
  }

  const int nsample = data_j.value(0).size();
  FitParameterDistribution params(nsample, guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);

  fit.fit(params, chisq, chisq_per_dof, data_j_inrange);

  double dof = chisq.sample(0)/chisq_per_dof.sample(0);
  jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });

  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "Dof: " << dof << std::endl;
  std::cout << "P-value: " << pvalue << std::endl;
  
#ifdef HAVE_HDF5
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(params, "params.hdf5"); 
  writeParamsStandard(pvalue, "pvalue.hdf5");
#endif

  FitFuncPolicy<FitFunc,ArgsType>::plot(args,*fitfunc,data_j,params);
  delete fitfunc;
}

#endif
