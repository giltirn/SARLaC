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

  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    int nsample = sigma(0).size();
    NumericSquareMatrix<jackknifeDistribution<double> > cov(sigma.size(), jackknifeDistribution<double>(nsample,0.));
    for(int i=0;i<sigma.size();i++) cov(i,i) = sigma(i) * sigma(i);      
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
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
    jackknifeDistribution<double> condition_number;
    svd_inverse(inv_corr, corr, condition_number);

    //Test the quality of the inverse
    NumericSquareMatrix<jackknifeDistribution<double> > test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistribution<double>(nsample,1.0);    
    jackknifeDistribution<double> resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number.mean() << " +- " << condition_number.standardError()/sqrt(nsample-1.) << std::endl;
    std::cout << "||CorrMat * CorrMat^{-1} - 1||_E = " << resid.mean() << " +- " << resid.standardError()/sqrt(nsample-1.) << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_corr,sigma);
  }

  void setUncorrelated(){ //because the fitter stores pointers we can modify the correlation matrix in place
    std::cout << "Setting correlation matrix to unit matrix\n";
    int nsample = sigma(0).size();
    for(int i=0;i<sigma.size();i++)
      for(int j=0;j<sigma.size();j++)
	inv_corr(i,j) = corr(i,j) = jackknifeDistribution<double>(nsample, i==j ? 1. : 0.);
  }

  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    NumericSquareMatrix<jackknifeDistribution<double> > cov = corr;
    for(int i=0;i<corr.size();i++)
      for(int j=0;j<corr.size();j++)
	cov(i,j) = cov(i,j) * sigma(i) * sigma(j);
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }
};



template<typename FitPolicies>
struct importCostFunctionParameters<correlatedCovFitPolicy,FitPolicies>{
  static_assert(std::is_same<typename FitPolicies::DistributionType,jackknifeDistribution<double> >::value, "Currently only support jackknifeDistribution<double>");
  
  NumericSquareMatrix<jackknifeDistribution<double> > cov;
  NumericSquareMatrix<jackknifeDistribution<double> > inv_cov;

  template<typename GeneralizedCoord, template<typename> class V>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, doubleJackknifeDistribution<double,V> > &data){

    const int nsample = data.value(0).size();
    const int ndata = data.size();   
    cov = NumericSquareMatrix<jackknifeDistribution<double> >(ndata);
    for(int i=0;i<ndata;i++){
      cov(i,i) = doubleJackknifeDistribution<double,V>::covariance(data.value(i), data.value(i));
      for(int j=i+1;j<ndata;j++)
	cov(i,j) = cov(j,i) = doubleJackknifeDistribution<double,V>::covariance(data.value(i),data.value(j));
    }
    
    inv_cov.resize(ndata, jackknifeDistribution<double>(nsample));
    jackknifeDistribution<double> condition_number;
    svd_inverse(inv_cov, cov, condition_number);

    //Test the quality of the inverse
    NumericSquareMatrix<jackknifeDistribution<double> > test = cov * inv_cov;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistribution<double>(nsample,1.0);    
    jackknifeDistribution<double> resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number.mean() << " +- " << condition_number.standardError()/sqrt(nsample-1.) << std::endl;
    std::cout << "||CovMat * CovMat^{-1} - 1||_E = " << resid.mean() << " +- " << resid.standardError()/sqrt(nsample-1.) << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_cov);
  }

  void setUncorrelated(){ //because the fitter stores pointers we can modify the correlation matrix in place
    std::cout << "Setting correlation matrix to unit matrix\n";
    int nsample = cov(0,0).size();
    for(int i=0;i<cov.size();i++){
      for(int j=0;j<cov.size();j++){
	if(i!=j){
	  inv_cov(i,j) = cov(i,j) = jackknifeDistribution<double>(nsample, 0.);
	}else{
	  inv_cov(i,j) = jackknifeDistribution<double>(nsample, 1.)/cov(i,j);
	}
      }
    }
  }

  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }
};



CPSFIT_END_NAMESPACE
#endif
