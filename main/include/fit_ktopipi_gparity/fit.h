#ifndef _KTOPIPI_FIT_H_
#define _KTOPIPI_FIT_H_

template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
std::vector<jackknifeDistribution<typename FitFunc::Params> > fit_corr_uncorr(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
									      const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_fit_dj,
									      const Args &args, const CMDline &cmdline){
  const int nsample = A0_fit_j[0].value(0).size();
  typedef typename FitFunc::Params Params;
  std::vector<Params> guess(10);
  std::vector<jackknifeDistribution<Params> > fit_params(10);

  //Perform the fit
  typedef typename composeFitPolicy<amplitudeDataCoord, FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
  FitFunc fitfunc;
  
  for(int q=0;q<10;q++){
    std::cout << "Starting fit for Q=" << q+1 << std::endl;
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);

    //Compute the correlation matrix / weights
    importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, A0_fit_dj[q]);

    readFrozenParams(fitter, q+1, cmdline, nsample);
  
    jackknifeDistribution<Params> &params = fit_params[q];
    params = jackknifeDistribution<Params>(nsample, guess[q]);
    jackknifeDistributionD chisq;
    jackknifeDistributionD chisq_per_dof;
    fitter.fit(params, chisq, chisq_per_dof, A0_fit_j[q]);

    distributionPrint<jackknifeDistribution<Params> >::printer(new ktopipiParamsPrinter<FitFunc>);

    std::cout << "Q" << q+1 << " Params: " << params << std::endl;
    std::cout << "Q" << q+1 << " Chisq: " << chisq << std::endl;
    std::cout << "Q" << q+1 << " Chisq/dof: " << chisq_per_dof << std::endl;
  }

  return fit_params;
}

template<typename FitFunc>
inline std::vector<jackknifeDistribution<typename FitFunc::Params> > fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
									 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_fit_dj,
									 const Args &args, const CMDline &cmdline){
  return args.correlated ?
    fit_corr_uncorr<FitFunc,correlatedFitPolicy>(A0_fit_j,A0_fit_dj,args,cmdline) :
    fit_corr_uncorr<FitFunc,uncorrelatedFitPolicy>(A0_fit_j,A0_fit_dj,args,cmdline);
}

#endif
