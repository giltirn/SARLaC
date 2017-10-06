#ifndef _KTOPIPI_FIT_H_
#define _KTOPIPI_FIT_H_

//Fit to each Q independently
template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
struct fit_corr_uncorr{
  static std::vector<jackknifeDistribution<typename FitFunc::Params> > fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
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
};

//Simultaneous fit
template<template<typename> class corrUncorrFitPolicy>
struct fit_corr_uncorr<FitKtoPiPiSim, corrUncorrFitPolicy>{
  static std::vector<jackknifeDistribution<FitKtoPiPiSim::Params> > fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
									const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_fit_dj,
									const Args &args, const CMDline &cmdline){
    const int nsample = A0_fit_j[0].value(0).size();
    typedef FitKtoPiPiSim FitFunc;
    typedef FitFunc::Params Params;
    Params guess;
    
    //Gather the data into single containers
    int sz = 0;
    for(int q=0;q<10;q++) sz += A0_fit_j[q].size();
    std::vector<std::pair<int,int> > map(sz);
    int a=0;
    for(int q=0;q<10;q++)
      for(int i=0;i<A0_fit_j[q].size();i++)
	map[a++] = std::pair<int,int>(q,i);
        
    correlationFunction<amplitudeDataCoordSim, jackknifeDistributionD> A0_fit_j_all(sz,
										    [&](const int s){
										      const amplitudeDataCoord &c = A0_fit_j[map[s].first].coord(map[s].second);
										      const jackknifeDistributionD &v = A0_fit_j[map[s].first].value(map[s].second);
										      amplitudeDataCoordSim cc(c, map[s].first);
										      return typename correlationFunction<amplitudeDataCoordSim, jackknifeDistributionD>::ElementType(cc,v);
										    });
    correlationFunction<amplitudeDataCoordSim, doubleJackknifeDistributionD> A0_fit_dj_all(sz,
											  [&](const int s){
											    const amplitudeDataCoord &c = A0_fit_dj[map[s].first].coord(map[s].second);
											    const doubleJackknifeDistributionD &v = A0_fit_dj[map[s].first].value(map[s].second);
											    amplitudeDataCoordSim cc(c, map[s].first);
											    return typename correlationFunction<amplitudeDataCoordSim, doubleJackknifeDistributionD>::ElementType(cc,v);
											  });
    //Perform the fit
    typedef typename composeFitPolicy<amplitudeDataCoordSim, FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
    FitFunc fitfunc;
  
    std::cout << "Starting fit for all Q simultaneously" << std::endl;
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);

    //Compute the correlation matrix / weights
    importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, A0_fit_dj_all);

    readFrozenParams(fitter, -1, cmdline, nsample);

    distributionPrint<jackknifeDistribution<Params> >::printer(new ktopipiParamsPrinter<FitFunc>);
    
    jackknifeDistribution<Params> params(nsample, guess);
    jackknifeDistributionD chisq;
    jackknifeDistributionD chisq_per_dof;
    fitter.fit(params, chisq, chisq_per_dof, A0_fit_j_all);

    std::cout << "Params: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

    return std::vector<jackknifeDistribution<typename FitFunc::Params> >(1, std::move(params));
  }
};





template<typename FitFunc>
inline std::vector<jackknifeDistribution<typename FitFunc::Params> > fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
									 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_fit_dj,
									 const Args &args, const CMDline &cmdline){
  return args.correlated ?
    fit_corr_uncorr<FitFunc,correlatedFitPolicy>::fit(A0_fit_j,A0_fit_dj,args,cmdline) :
    fit_corr_uncorr<FitFunc,uncorrelatedFitPolicy>::fit(A0_fit_j,A0_fit_dj,args,cmdline);
}

void getFitData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
		std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_fit_dj,		
		const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_j,
		const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_dj,
		const Args &args){
  //Separate out the data in the desired fit range
  for(int q=0;q<10;q++){
    for(int d=0;d<A0_j[q].size();d++){
      const int t = int(A0_j[q].coord(d).t);
      const int tsep_k_pi = A0_j[q].coord(d).tsep_k_pi;
      const int tsep_op_pi = tsep_k_pi - t;
      if(t <= tsep_k_pi && t >= args.tmin_k_op && tsep_op_pi >= args.tmin_op_pi){
	A0_fit_j[q].push_back(A0_j[q].coord(d), A0_j[q].value(d));
	A0_fit_dj[q].push_back(A0_dj[q].coord(d), A0_dj[q].value(d));
      }
    }
  }
}


template<typename FitFunc>
inline void fitAndPlotFF(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
			 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_all_dj,
			 const Args &args, const CMDline &cmdline){
  //Extract the data we are going to fit
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_fit_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > A0_fit_dj(10);
  getFitData(A0_fit_j,A0_fit_dj,A0_all_j,A0_all_dj,args);
  
  std::cout << "Including " << A0_fit_j[0].size() << " data points in fit\n";
  for(int q=0;q<10;q++){
    std::cout << "For Q" << q+1 << std::endl;
    for(int i=0;i<A0_fit_j[q].size();i++)
      std::cout << A0_fit_j[q].coord(i) << " : " << A0_fit_j[q].value(i) << std::endl;
  }
  
  std::vector<jackknifeDistribution<typename FitFunc::Params> > fit_params = fit<FitFunc>(A0_fit_j,A0_fit_dj,args,cmdline);

  extractMdata<FitFunc> extractor(fit_params);
  plotErrorWeightedData(A0_all_j,extractor,args);

#ifdef HAVE_HDF5
  writeParamsStandard(fit_params, "params.hdf5");
#endif
}

inline void fitAndPlot(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
		       const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_all_dj,
		       const Args &args, const CMDline &cmdline){
  switch(args.fitfunc){
  case FitSeparate:
    return fitAndPlotFF<FitKtoPiPi>(A0_all_j,A0_all_dj,args,cmdline);
  case FitSimultaneous:
    return fitAndPlotFF<FitKtoPiPiSim>(A0_all_j,A0_all_dj,args,cmdline);
  default:
    error_exit(std::cout << "fitAndPlot(..) Unknown fit function " << args.fitfunc << std::endl);
  }
}
  

#endif
