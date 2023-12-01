#ifndef FIT_SIMULTANEOUS_FIT_H
#define FIT_SIMULTANEOUS_FIT_H

typedef FitFuncSimultaneous FitFunc;

template<template<typename> class CostFunctionPolicy>
void fitSpecCorr(const jackknifeSimFitCorrelationFunction &data_j, const doubleJackknifeSimFitCorrelationFunction data_dj, const FitFunc &fitfunc, const Args &args){
  jackknifeSimFitCorrelationFunction data_j_inrange;
  doubleJackknifeSimFitCorrelationFunction data_dj_inrange;
  for(int i=0;i<data_j.size();i++){
    int t = (int)data_j.coord(i).t;
    if(t >= args.t_min && t <= args.t_max){
      data_j_inrange.push_back(data_j[i]);
      data_dj_inrange.push_back(data_dj[i]);
    }
  }

  typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, CostFunctionPolicy>::type FitPolicies;
  typedef fitter<FitPolicies> Fitter;
  typedef typename FitPolicies::FitParameterDistribution FitParameterDistribution;
  
  Fitter fit;
  fit.importFitFunc(fitfunc);

  importCostFunctionParameters<CostFunctionPolicy,FitPolicies> prepare(fit,data_dj_inrange);

  typename FitFunc::ParameterType guess = fitfunc.guess();
  // if(cmdline.load_guess){
  //   parse(guess,cmdline.guess_file);
  // }

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
}


void fit(const jackknifeSimFitCorrelationFunction &data_j, const doubleJackknifeSimFitCorrelationFunction data_dj, const FitFunc &fitfunc, const Args &args){
  return args.correlated ? 
    fitSpecCorr<correlatedFitPolicy>(data_j,data_dj,fitfunc,args) : 
    fitSpecCorr<uncorrelatedFitPolicy>(data_j, data_dj,fitfunc,args);
}


#endif
