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
  typedef typename costFunctionType::ValueType ValueType;

  //To ensure each loop iteration independent we can't store global state, instead we have thread local state
  class costFunctionState: public fitFuncPolicyState{
  private:
    sampleSeriesType data_s;
    sampleInvCorrType inv_corr_s;
    std::vector<double> sigma_s;
    std::unique_ptr<costFunctionType> cost_func;

  public:
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const CorrelationFunctionDistribution &data, const MatrixDistribution &inv_corr, const VectorDistribution &sigma, 
		      const std::vector<typename costFunctionType::Prior> &priors, const int s):
      data_s(data,s), inv_corr_s(inv_corr,s), sigma_s(data.size()), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      for(int d=0;d<data.size();d++)
	sigma_s[d] = iterate<DistributionType>::at(s,sigma[d]);
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma_s, inv_corr_s));
      for(int p=0;p<priors.size();p++) cost_func->addPrior(priors[p]);
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    int Ndof() const{ return cost_func->Ndof(); }
  };

private:
  MatrixDistribution const* inv_corr;
  VectorDistribution const* sigma;
  
  std::vector<typename costFunctionType::Prior> priors;

public:
  correlatedFitPolicy(): inv_corr(NULL), sigma(NULL){}
  
  //Import the inverse correlation matrix and sigma
  void importCostFunctionParameters(const MatrixDistribution &_inv_corr, 
				    const VectorDistribution &_sigma){
    inv_corr = &_inv_corr;
    sigma = &_sigma;
  }

  //Add a prior for constrained curve fitting (https://arxiv.org/pdf/hep-lat/0110175.pdf). This adds a term to chi^2 for the parameter p with index 'param_idx'
  //of the form d\Chi^2 = ( p - value )^2/weight^2
  void addPrior(ValueType value, ValueType weight, int param_idx){
    priors.push_back(typename costFunctionType::Prior(value,weight,param_idx));
  }

protected:
  inline costFunctionState generateCostFunctionState(const CorrelationFunctionDistribution &data, const int s){
    assert(inv_corr != NULL && sigma != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*inv_corr,*sigma,priors,s);
  }
};


template<typename FitFuncPolicy>
class correlatedCovFitPolicy: public FitFuncPolicy{
public:
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FitFuncPolicy);
  typedef CorrelatedCovChisqCostFunction<fitFunc, sampleSeriesType, sampleInvCorrType> costFunctionType;
  
  //To ensure each loop iteration independent we can't store global state, instead we have thread local state
  class costFunctionState: public fitFuncPolicyState{
  private:
    sampleSeriesType data_s;
    sampleInvCorrType inv_cov_s;
    std::unique_ptr<costFunctionType> cost_func;

  public:
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const CorrelationFunctionDistribution &data, const MatrixDistribution &inv_cov, const int s):
    data_s(data,s), inv_cov_s(inv_cov,s), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, inv_cov_s));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    int Ndof() const{ return data_s.size() - this->getFitFunc().Nparams(); }
  };

private:
  MatrixDistribution const* inv_cov;
public:
  correlatedCovFitPolicy(): inv_cov(NULL){}
  
  //Import the inverse correlation matrix and sigma
  void importCostFunctionParameters(const MatrixDistribution &_inv_cov){
    inv_cov = &_inv_cov;
  }
protected:
  inline costFunctionState generateCostFunctionState(const CorrelationFunctionDistribution &data, const int s){
    assert(inv_cov != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*inv_cov,s);
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
