#ifndef _SARLAC_FIT_WRAPPER_COSTFUNC_CORRELATED_COV_FIT_POLICY_H_
#define _SARLAC_FIT_WRAPPER_COSTFUNC_CORRELATED_COV_FIT_POLICY_H_

#include<config.h>
#include<utils/macros.h>
#include<fit/cost_function.h>
#include<fit/fit_wrapper/base_typedefs.h>
#include<fit/fit_wrapper/fitfunc_policy.h>

SARLAC_START_NAMESPACE

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

SARLAC_END_NAMESPACE
#endif
