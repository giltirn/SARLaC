#ifndef _CPSFIT_FIT_WRAPPER_COSTFUNC_IMPORT_H_
#define _CPSFIT_FIT_WRAPPER_COSTFUNC_IMPORT_H_

//The generation of the inverse correlation matrix (if applicable) and sigma for the cost function using the double-jackknife is something that often crops up. Here we implement a convenience wrapper to perform this boilerplate

#include<config.h>
#include<utils/macros.h>
#include<distribution/double_jackknife.h>
#include<fit/fit_wrapper/costfunction_policy.h>
#include<fit/fit_wrapper/fitter.h>

CPSFIT_START_NAMESPACE

template<template<typename> class corrUncorrFitPolicy, typename FitPolicies>
struct importCostFunctionParameters{};

//For uncorrelated fits
template<typename FitPolicies>
struct importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies>{
  static_assert(std::is_same<typename FitPolicies::DistributionType,jackknifeDistribution<double> >::value, "Currently only support jackknifeDistribution<double>");
  
  NumericVector<jackknifeDistribution<double> > sigma;

  template<typename GeneralizedCoord, template<typename> class V>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, doubleJackknifeDistribution<double,V> > &data): sigma(data.size()){
    for(int d=0;d<data.size();d++)
      sigma(d) = jackknifeDistribution<double>(sqrt(doubleJackknifeDistribution<double,V>::covariance(data.value(d) , data.value(d) ) ) );
    
    fitter.importCostFunctionParameters(sigma);
  }
};

//For correlated fits
template<typename FitPolicies>
struct importCostFunctionParameters<correlatedFitPolicy,FitPolicies>{
  static_assert(std::is_same<typename FitPolicies::DistributionType,jackknifeDistribution<double> >::value, "Currently only support jackknifeDistribution<double>");
  
  NumericSquareMatrix<jackknifeDistribution<double> > corr;
  NumericSquareMatrix<jackknifeDistribution<double> > inv_corr;
  NumericVector<jackknifeDistribution<double> > sigma;

  template<typename GeneralizedCoord, template<typename> class V>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, doubleJackknifeDistribution<double,V> > &data): sigma(data.size()){

    const int nsample = data.value(0).size();
    const int ndata = data.size();    
    NumericSquareMatrix<jackknifeDistribution<double>> cov(ndata);
    for(int i=0;i<ndata;i++){
      cov(i,i) = doubleJackknifeDistribution<double,V>::covariance(data.value(i), data.value(i));
      sigma(i) = sqrt(cov(i,i));

      for(int j=i+1;j<ndata;j++)
	cov(i,j) = cov(j,i) = doubleJackknifeDistribution<double,V>::covariance(data.value(i),data.value(j));
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


CPSFIT_END_NAMESPACE
#endif
