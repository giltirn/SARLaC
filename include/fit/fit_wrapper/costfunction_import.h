#ifndef _CPSFIT_FIT_WRAPPER_COSTFUNC_IMPORT_H_
#define _CPSFIT_FIT_WRAPPER_COSTFUNC_IMPORT_H_

//The generation of the inverse correlation matrix (if applicable) and sigma for the cost function using the double-jackknife is something that often crops up. Here we implement a convenience wrapper to perform this boilerplate

#include<config.h>
#include<utils/macros.h>
#include<distribution/double_jackknife.h>
#include<distribution/block_double_jackknife.h>
#include<fit/fit_wrapper/costfunction_policy.h>
#include<fit/fit_wrapper/fitter.h>

CPSFIT_START_NAMESPACE

template<template<typename> class corrUncorrFitPolicy, typename FitPolicies>
struct importCostFunctionParameters{};

//For uncorrelated fits
template<typename FitPolicies>
struct importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies>{
  typedef typename FitPolicies::DistributionType DistributionType;
  
  NumericVector<DistributionType> sigma;

  template<typename GeneralizedCoord, typename CovDistributionType> //CovDistributionType should have covariance static member
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, CovDistributionType> &data): sigma(data.size()){
    for(int d=0;d<data.size();d++)
      sigma(d) = sqrt(CovDistributionType::covariance(data.value(d) , data.value(d) ) );
    
    fitter.importCostFunctionParameters(sigma);
  }
  template<typename GeneralizedCoord, typename CovDistributionType>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, DistributionType > &data_j,
			       const correlationFunction<GeneralizedCoord, CovDistributionType > &data_dj):
    importCostFunctionParameters(fitter, data_dj){}

  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    NumericSquareMatrix<DistributionType> cov(sigma.size());					      
    for(int i=0;i<sigma.size();i++) cov(i,i) = sigma(i) * sigma(i);      
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }

};

//For correlated fits
template<typename FitPolicies>
struct importCostFunctionParameters<correlatedFitPolicy,FitPolicies>{
  typedef typename FitPolicies::DistributionType DistributionType;
  
  NumericSquareMatrix<DistributionType> corr;
  NumericSquareMatrix<DistributionType> inv_corr;
  NumericVector<DistributionType> sigma;

  template<typename GeneralizedCoord, typename CovDistributionType>
  void import(fitter<FitPolicies> &fitter,
	      const correlationFunction<GeneralizedCoord, CovDistributionType> &data){
    const int ndata = data.size();    
    sigma.resize(ndata);
    NumericSquareMatrix<DistributionType> cov(ndata);
    for(int i=0;i<ndata;i++){
      cov(i,i) = CovDistributionType::covariance(data.value(i), data.value(i));
      sigma(i) = sqrt(cov(i,i));

      for(int j=i+1;j<ndata;j++)
	cov(i,j) = cov(j,i) = CovDistributionType::covariance(data.value(i),data.value(j));
    }

    auto init = sigma(0).getInitializer();
    DistributionType one(init,1.);

    corr =  NumericSquareMatrix<DistributionType>(ndata);
    
    for(int i=0;i<ndata;i++){
      corr(i,i) = one;
      
      for(int j=i+1;j<ndata;j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
    }
    
    inv_corr = corr;
    DistributionType condition_number;
    svd_inverse(inv_corr, corr, condition_number);

    //Test the quality of the inverse
    NumericSquareMatrix<DistributionType> test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - one;
    DistributionType resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number.mean() << " +- " << condition_number.standardDeviation() << std::endl;
    std::cout << "||CorrMat * CorrMat^{-1} - 1||_E = " << resid.mean() << " +- " << resid.standardDeviation() << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_corr,sigma);
  }

  importCostFunctionParameters(){}

  template<typename GeneralizedCoord, typename CovDistributionType>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, CovDistributionType> &data){
    import(fitter, data);
  }

  template<typename GeneralizedCoord, template<typename> class V, typename CovDistributionType>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, DistributionType > &data_j,
			       const correlationFunction<GeneralizedCoord, CovDistributionType> &data_dj):
    importCostFunctionParameters(fitter, data_dj){}

  void setUncorrelated(){ //because the fitter stores pointers we can modify the correlation matrix in place
    std::cout << "Setting correlation matrix to unit matrix\n";
    DistributionType one(sigma(0).getInitializer(),1.);
    DistributionType zero(sigma(0).getInitializer(),0.);

    for(int i=0;i<sigma.size();i++)
      for(int j=0;j<sigma.size();j++)
	inv_corr(i,j) = corr(i,j) = i==j ? one : zero;
  }

  //Apply the Ledoit-Wolf procedure to stabilize the correlation matrix
  //(https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/chulwoo/notes/LedoitWolf/LedoitWolf.pdf)
  //Returns lambda
  double LedoitWolfShrinkage(){
    double lambda_num = 0.;
    double lambda_den = 0.;
    for(int i=0;i<corr.size();i++){
      for(int j=0;j<corr.size();j++){
	if(i != j){ double s = corr(i,j).standardError(); lambda_num += s*s; } //var(C(i,j))
	DistributionType c2 = corr(i,j)*corr(i,j);
	lambda_den += c2.mean();
      }
    }
    double lambda = lambda_num/lambda_den;

    std::cout << "Applying Ledoit-Wolf procedure with lambda = " << lambda << std::endl;
    
    auto init = sigma(0).getInitializer();

    DistributionType lambda_j(init, lambda);
    DistributionType _1mlambda_j(init, 1. - lambda);

    for(int i=0;i<corr.size();i++){
      for(int j=0;j<corr.size();j++){
	if(i==j) corr(i,j) = lambda_j + _1mlambda_j * corr(i,j);
	else corr(i,j) = _1mlambda_j * corr(i,j);
      }
    }
    DistributionType condition_number;
    svd_inverse(inv_corr, corr, condition_number);

    //Test the quality of the inverse
    DistributionType one(init, 1.);
    NumericSquareMatrix<DistributionType> test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - one;    
    DistributionType resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number.mean() << " +- " << condition_number.standardDeviation() << std::endl;
    std::cout << "||CorrMat * CorrMat^{-1} - 1||_E = " << resid.mean() << " +- " << resid.standardDeviation() << std::endl;
    return lambda;
  }



  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    NumericSquareMatrix<DistributionType> cov = corr;
    for(int i=0;i<corr.size();i++)
      for(int j=0;j<corr.size();j++)
	cov(i,j) = cov(i,j) * sigma(i) * sigma(j);
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }
};

//Note    Corr(i,j) = Cov(i,j)/\sigma(i)/\sigma(j)   thus   Corr = \Sigma^{-1} Cov \Sigma^{-1}   where \Sigma = diag(\sigma(i))
//Hence  Cov = \Sigma Corr \Sigma   and  Cov^{-1} = \Sigma^{-1} Corr^{-1} \Sigma^{-1}
//Removing \Sigma conditions the Correlation matrix so it is always wise to do

template<typename FitPolicies>
struct importCostFunctionParameters<correlatedCovFitPolicy,FitPolicies>{
  typedef typename FitPolicies::DistributionType DistributionType;
  
  NumericSquareMatrix<DistributionType> cov;
  NumericSquareMatrix<DistributionType> inv_cov;

  importCostFunctionParameters<correlatedFitPolicy,FitPolicies> covp;

  importCostFunctionParameters(): covp(NULL){}

  template<typename GeneralizedCoord, typename CovDistributionType>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, CovDistributionType> &data){
    int ndata = data.size();
    covp.importCostFunctionParameters(fitter, data);
    
    cov = covp.corr;
    inv_cov = covp.inv_corr;
    for(int i=0;i<ndata;i++)
      for(int j=0;j<ndata;j++){
	cov(i,j) = cov(i,j) * covp.sigma(i) * covp.sigma(j);
	inv_cov(i,j) = inv_cov(i,j) / covp.sigma(i) / covp.sigma(j);
      }
  }

  template<typename GeneralizedCoord, typename CovDistributionType>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, DistributionType > &data_j,
			       const correlationFunction<GeneralizedCoord, CovDistributionType > &data_dj):
    importCostFunctionParameters(fitter, data_dj){}

  void setUncorrelated(){ //because the fitter stores pointers we can modify the correlation matrix in place
    std::cout << "Setting correlation matrix to unit matrix\n";
    auto init = cov(0,0).getInitializer();
    DistributionType one(init,1.);
    DistributionType zero(init,0.);

    for(int i=0;i<cov.size();i++){
      for(int j=0;j<cov.size();j++){
	if(i!=j){
	  inv_cov(i,j) = cov(i,j) = zero;
	}else{
	  inv_cov(i,j) = one /cov(i,j);
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


//For frozen correlated fits
template<typename FitPolicies>
struct importCostFunctionParameters<frozenCorrelatedFitPolicy,FitPolicies>{
  typedef typename FitPolicies::DistributionType DistributionType;
  
  NumericSquareMatrix<double> corr;
  NumericSquareMatrix<double> inv_corr;
  std::vector<double> sigma;

  template<typename GeneralizedCoord>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, DistributionType> &data): sigma(data.size()){

    const int ndata = data.size();    
    NumericSquareMatrix<double> cov(ndata);
    for(int i=0;i<ndata;i++){
      cov(i,i) = DistributionType::covariance(data.value(i), data.value(i));
      sigma[i] = sqrt(cov(i,i));

      for(int j=i+1;j<ndata;j++)
	cov(i,j) = cov(j,i) = DistributionType::covariance(data.value(i),data.value(j));
    }

    corr =  NumericSquareMatrix<double>(ndata);
    
    for(int i=0;i<ndata;i++){
      corr(i,i) = 1.;
      
      for(int j=i+1;j<ndata;j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma[i]/sigma[j];
    }
    
    inv_corr = corr;
    double condition_number;
    svd_inverse(inv_corr, corr, condition_number);

    //Test the quality of the inverse
    NumericSquareMatrix<double> test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - 1.0;    
    double resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number << std::endl;
    std::cout << "||CorrMat * CorrMat^{-1} - 1||_E = " << resid << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_corr,sigma);
  }
  
  //To keep the interface the same as the others, we allow for the passing of a second correlation function containing data that would be used to obtain the unfrozen covariance matrix, but it is not used here
  template<typename GeneralizedCoord, typename UnusedCovDistributionType>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, DistributionType> &data_j,
			       const correlationFunction<GeneralizedCoord, UnusedCovDistributionType> &data_dj):
    importCostFunctionParameters(fitter, data_j){}

  void setUncorrelated(){ //because the fitter stores pointers we can modify the correlation matrix in place
    std::cout << "Setting correlation matrix to unit matrix\n";
    for(int i=0;i<sigma.size();i++)
      for(int j=0;j<sigma.size();j++)
	inv_corr(i,j) = corr(i,j) = i==j ? 1. : 0.;
  }

  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    NumericSquareMatrix<double> cov = corr;
    for(int i=0;i<corr.size();i++)
      for(int j=0;j<corr.size();j++)
	cov(i,j) = cov(i,j) * sigma[i] * sigma[j];
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }
};



CPSFIT_END_NAMESPACE
#endif
