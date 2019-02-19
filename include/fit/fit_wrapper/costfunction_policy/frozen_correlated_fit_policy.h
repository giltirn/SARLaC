#ifndef _CPSFIT_FIT_WRAPPER_COSTFUNC_FROZEN_CORRELATED_FIT_POLICY_H_
#define _CPSFIT_FIT_WRAPPER_COSTFUNC_FROZEN_CORRELATED_FIT_POLICY_H_

#include<config.h>
#include<utils/macros.h>
#include<fit/cost_function.h>
#include<fit/fit_wrapper/base_typedefs.h>
#include<fit/fit_wrapper/fitfunc_policy.h>

#include "macros.h"

CPSFIT_START_NAMESPACE

//Use a frozen correlation matrix (i.e. same matrix for each sample) obtained from the data jackknife
template<typename FitFuncPolicy>
class frozenCorrelatedFitPolicy: public FitFuncPolicy{
public:
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FitFuncPolicy);
  typedef NumericSquareMatrix<double> InvCorrType;
  typedef std::vector<double> SigmaType;

  typedef CorrelatedChisqCostFunction<fitFunc, sampleSeriesType, InvCorrType> costFunctionType;
  typedef typename costFunctionType::ValueType ValueType;

  //To ensure each loop iteration independent we can't store global state, instead we have thread local state
  class costFunctionState: public fitFuncPolicyState{
  private:
    sampleSeriesType data_s;
    std::unique_ptr<costFunctionType> cost_func;

  public:
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const CorrelationFunctionDistribution &data, const InvCorrType &inv_corr, const SigmaType &sigma, const int s):
      data_s(data,s), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma, inv_corr));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    int Ndof() const{ return cost_func->Ndof(); }
  };

private:
  InvCorrType const* inv_corr;
  SigmaType const* sigma;

public:
  frozenCorrelatedFitPolicy(): inv_corr(NULL), sigma(NULL){}
  
  //Import the inverse correlation matrix and sigma
  void importCostFunctionParameters(const InvCorrType &_inv_corr, 
				    const SigmaType &_sigma){
    inv_corr = &_inv_corr;
    sigma = &_sigma;
  }

protected:
  inline costFunctionState generateCostFunctionState(const CorrelationFunctionDistribution &data, const int s){
    assert(inv_corr != NULL && sigma != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*inv_corr,*sigma,s);
  }
};

CPSFIT_END_NAMESPACE
#endif
