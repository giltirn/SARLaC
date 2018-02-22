#ifndef _CPSFIT_FIT_WRAPPER_COSTFUNC_POLICY_H_
#define _CPSFIT_FIT_WRAPPER_COSTFUNC_POLICY_H_

//The cost function policy is responsible for handling the cost function used in the minimizer. It inherits methods from the fit function policy.

#include<config.h>
#include<utils/macros.h>
#include<fit/cost_function.h>
#include<fit/fit_wrapper/base_typedefs.h>
#include<fit/fit_wrapper/fitfunc_policy.h>


CPSFIT_START_NAMESPACE

#define INHERIT_FIT_POLICY_TYPEDEFS(FROM) \
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FROM); \
  INHERIT_TYPEDEF(FROM,costFunctionType); \
  INHERIT_TYPEDEF(FROM,costFunctionState)


template<typename FitFuncPolicy>
class correlatedFitPolicy: public FitFuncPolicy{
public:
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FitFuncPolicy);
  typedef CorrelatedChisqCostFunction<fitFunc, sampleSeriesType, sampleInvCorrType> costFunctionType;
  
  //To ensure each loop iteration independent we can't store global state, instead we have thread local state
  class costFunctionState: public fitFuncPolicyState{
  private:
    sampleSeriesType data_s;
    sampleInvCorrType inv_corr_s;
    std::vector<double> sigma_s;
    std::unique_ptr<costFunctionType> cost_func;

  public:
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const CorrelationFunctionDistribution &data, const MatrixDistribution &inv_corr, const VectorDistribution &sigma, const int s):
      data_s(data,s), inv_corr_s(inv_corr,s), sigma_s(data.size()), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      for(int d=0;d<data.size();d++)
	sigma_s[d] = iterate<DistributionType>::at(s,sigma[d]);
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma_s, inv_corr_s));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    int Ndof() const{ return data_s.size() - this->getFitFunc().Nparams(); }
  };

private:
  MatrixDistribution const* inv_corr;
  VectorDistribution const* sigma;
public:
  correlatedFitPolicy(): inv_corr(NULL), sigma(NULL){}
  
  //Import the inverse correlation matrix and sigma
  void importCostFunctionParameters(const MatrixDistribution &_inv_corr, 
				    const VectorDistribution &_sigma){
    inv_corr = &_inv_corr;
    sigma = &_sigma;
  }
protected:
  inline costFunctionState generateCostFunctionState(const CorrelationFunctionDistribution &data, const int s){
    assert(inv_corr != NULL && sigma != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*inv_corr,*sigma,s);
  }
};



template<typename FitFuncPolicy>
class uncorrelatedFitPolicy: public FitFuncPolicy{
public:
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FitFuncPolicy);
  typedef UncorrelatedChisqCostFunction<fitFunc, sampleSeriesType> costFunctionType;
  
  //To ensure each loop iteration independent we can't store global state, instead we have thread local state
  class costFunctionState: public fitFuncPolicyState{
  private:
    sampleSeriesType data_s;
    std::vector<double> sigma_s;
    std::unique_ptr<costFunctionType> cost_func;

  public:
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const CorrelationFunctionDistribution &data, const VectorDistribution &sigma, const int s):
      data_s(data,s), sigma_s(data.size()), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      for(int d=0;d<data.size();d++)
	sigma_s[d] = iterate<DistributionType>::at(s,sigma[d]);
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma_s));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    int Ndof() const{ return data_s.size() - this->getFitFunc().Nparams(); }
  };

private:
  VectorDistribution const* sigma;
public:
  uncorrelatedFitPolicy(): sigma(NULL){}
  
  void importCostFunctionParameters(const VectorDistribution &_sigma){
    sigma = &_sigma;
  }
protected:
  inline costFunctionState generateCostFunctionState(const CorrelationFunctionDistribution &data, const int s){
    assert(sigma != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*sigma,s);
  }
};



CPSFIT_END_NAMESPACE
#endif
