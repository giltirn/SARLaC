#ifndef _KTOPIPI_FIT_SAMA_EXPAND_H_
#define  _KTOPIPI_FIT_SAMA_EXPAND_H_

//Expand the SAMA jackknife by appending NQ samples
jackknifeDistributionD SAMAexpandJack(const jackknifeDistributionD &in, const int NQ){
  double mean = in.mean();
  int N = in.size();
  jackknifeDistributionD out(N+NQ);
  for(int i=0;i<N;i++) out.sample(i) = in.sample(i);
  for(int i=0;i<NQ;i++) out.sample(N+i) = mean;	
  return out;
}


template<template<typename> class corrUncorrFitPolicy, typename FitPolicies>
struct importCostFunctionParametersSAMAexpand{};

//For uncorrelated fits
template<typename FitPolicies>
struct importCostFunctionParametersSAMAexpand<uncorrelatedFitPolicy,FitPolicies>{
  static_assert(std::is_same<typename FitPolicies::DistributionType,jackknifeDistribution<double> >::value, "Currently only support jackknifeDistribution<double>");
  
  NumericVector<jackknifeDistribution<double> > sigma;

  template<typename GeneralizedCoord, template<typename> class V>
  importCostFunctionParametersSAMAexpand(fitter<FitPolicies> &fitter,
					 const correlationFunction<GeneralizedCoord, doubleJackknifeDistribution<double,V> > &data, const int NQ): sigma(data.size()){
    for(int d=0;d<data.size();d++)
      sigma(d) = SAMAexpandJack(jackknifeDistribution<double>(sqrt(doubleJackknifeDistribution<double,V>::covariance(data.value(d) , data.value(d) ) ) ), NQ);
    
    fitter.importCostFunctionParameters(sigma);
  }
};

//For correlated fits
template<typename FitPolicies>
struct importCostFunctionParametersSAMAexpand<correlatedFitPolicy,FitPolicies>{
  static_assert(std::is_same<typename FitPolicies::DistributionType,jackknifeDistribution<double> >::value, "Currently only support jackknifeDistribution<double>");
  
  NumericSquareMatrix<jackknifeDistribution<double> > corr;
  NumericSquareMatrix<jackknifeDistribution<double> > inv_corr;
  NumericVector<jackknifeDistribution<double> > sigma;

  template<typename GeneralizedCoord, template<typename> class V>
  importCostFunctionParametersSAMAexpand(fitter<FitPolicies> &fitter,
					 const correlationFunction<GeneralizedCoord, doubleJackknifeDistribution<double,V> > &data, const int NQ): sigma(data.size()){

    const int nsample = data.value(0).size() + NQ;
    const int ndata = data.size();    
    NumericSquareMatrix<jackknifeDistribution<double>> cov(ndata);
    for(int i=0;i<ndata;i++){
      cov(i,i) = SAMAexpandJack(doubleJackknifeDistribution<double,V>::covariance(data.value(i), data.value(i)), NQ);
      sigma(i) = sqrt(cov(i,i));

      for(int j=i+1;j<ndata;j++)
	cov(i,j) = cov(j,i) = SAMAexpandJack(doubleJackknifeDistribution<double,V>::covariance(data.value(i),data.value(j)),NQ);
    }

    corr =  NumericSquareMatrix<jackknifeDistribution<double> >(ndata);
    
    for(int i=0;i<ndata;i++){
      corr(i,i) = jackknifeDistribution<double>(nsample,1.);
      
      for(int j=i+1;j<ndata;j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
    }
    
    inv_corr.resize(ndata, jackknifeDistribution<double>(nsample));
    svd_inverse(inv_corr, corr);

    //Test the quality of the inverse
    NumericSquareMatrix<jackknifeDistribution<double> > test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistribution<double>(nsample,1.0);    
    std::cout << "|CorrMat * CorrMat^{-1} - 1|^2 = " << mod2(test) << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_corr,sigma);
  }
};



//Fit to each Q independently
template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
struct fit_corr_uncorr_SAMAexpand{
  static std::vector<jackknifeDistribution<typename FitFunc::Params> > fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
									   const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,
									   const SampleAMAargs &args, const SampleAMAcmdLine &cmdline){
    const int nsample = A0_fit_j[0].value(0).size(); //SAMA expanded already
    typedef typename FitFunc::Params Params;
    std::vector<Params> guess(10);
    std::vector<jackknifeDistribution<Params> > fit_params(10);

    //Perform the fit
    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
    FitFunc fitfunc;
  
    for(int q=0;q<10;q++){
      std::cout << "Starting fit for Q=" << q+1 << std::endl;
      fitter<FitPolicies> fitter;
      fitter.importFitFunc(fitfunc);

      //Compute the correlation matrix / weights
      importCostFunctionParametersSAMAexpand<corrUncorrFitPolicy,FitPolicies> prepare(fitter, A0_fit_dj[q], cmdline.SAMAexpandN);

      if(cmdline.load_freeze_data)
	readFrozenParams(fitter,q+1,cmdline.freeze_data,nsample);
  
      jackknifeDistribution<Params> &params = fit_params[q];
      params = jackknifeDistribution<Params>(nsample, guess[q]);
      jackknifeDistributionD chisq;
      jackknifeDistributionD chisq_per_dof;
      fitter.fit(params, chisq, chisq_per_dof, A0_fit_j[q]);

      int nparams = params.sample(0).size();
      for(int p=0;p<nparams;p++){
	jackknifeDistributionD tmp(nsample, [&](const int s){ return params.sample(s)(p); }); 
	std::cout << params.sample(0).memberName(p) << ": " << tmp << std::endl;
      }

      std::cout << "Q" << q+1 << " Chisq: " << chisq << std::endl;
      std::cout << "Q" << q+1 << " Chisq/dof: " << chisq_per_dof << std::endl;
    }

    return fit_params;
  }
};



template<typename FitFunc>
inline typename fitReturnType<FitFunc>::type fitSAMAexpand(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
						 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,
						 const SampleAMAargs &args, const SampleAMAcmdLine &cmdline){
  return args.correlated ?
    fit_corr_uncorr_SAMAexpand<FitFunc,correlatedFitPolicy>::fit(A0_fit_j,A0_fit_dj,args,cmdline) :
    fit_corr_uncorr_SAMAexpand<FitFunc,uncorrelatedFitPolicy>::fit(A0_fit_j,A0_fit_dj,args,cmdline);
}




template<typename FitFunc>
inline void fitAndPlotFFSAMAexpand(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
			 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
			 const SampleAMAargs &args, const SampleAMAcmdLine &cmdline){

  //SAMA expand data
  int NQ = cmdline.SAMAexpandN;
  typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD>::ElementType Elem;
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j_exp(A0_all_j.size());
  for(int i=0;i<A0_all_j.size();i++)
    A0_all_j_exp[i].resize(A0_all_j[i].size(),
			   [&](const int e){ return Elem(A0_all_j[i][e].first, SAMAexpandJack(A0_all_j[i][e].second,NQ)); }
			   );
  
  if(A0_all_j_exp[0].value(0).size() != A0_all_j[0].value(0).size() + NQ)
    error_exit(std::cout << "SAMA failed, got " << A0_all_j_exp[0].value(0).size() << " expected " << A0_all_j[0].value(0).size() + NQ << std::endl);

  //Extract the data we are going to fit
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_fit_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_fit_dj(10);
  getFitData(A0_fit_j,A0_fit_dj,A0_all_j_exp,A0_all_dj, args.tmin_k_op, args.tmin_op_pi);
  
  std::cout << "Including " << A0_fit_j[0].size() << " data points in fit\n";
  for(int q=0;q<10;q++){
    std::cout << "For Q" << q+1 << std::endl;
    for(int i=0;i<A0_fit_j[q].size();i++)
      std::cout << A0_fit_j[q].coord(i) << " : " << A0_fit_j[q].value(i) << std::endl;
  }
  
  typename fitReturnType<FitFunc>::type fit_params = fitSAMAexpand<FitFunc>(A0_fit_j,A0_fit_dj,args,cmdline);

  extractMdata<FitFunc> extractor(fit_params);
  plotErrorWeightedData(A0_all_j_exp,extractor,args.tmin_k_op,args.tmin_op_pi);

#ifdef HAVE_HDF5
  writeParamsStandard(fit_params, "params.hdf5");
#endif
}

inline void fitAndPlotSAMAexpand(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
				 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
				 const SampleAMAargs &args, const SampleAMAcmdLine &cmdline){
  switch(args.fitfunc){
  case KtoPiPiFitFunc::FitSeparate:
    return fitAndPlotFFSAMAexpand<FitKtoPiPi>(A0_all_j,A0_all_dj,args,cmdline);
  case KtoPiPiFitFunc::FitSeparateTwoExp:
    return fitAndPlotFFSAMAexpand<FitKtoPiPiTwoExp>(A0_all_j,A0_all_dj,args,cmdline);
  default:
    error_exit(std::cout << "fitAndPlotSAMAexpand(..) Unknown fit function " << args.fitfunc << std::endl);
  }
}



#endif
