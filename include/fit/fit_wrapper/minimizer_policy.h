#ifndef _CPSFIT_FIT_WRAPPER_MINIMIZER_POLICY_H
#define _CPSFIT_FIT_WRAPPER_MINIMIZER_POLICY_H

//The minimizer policy is the topmost layer, and defines which algorithm is used to minimize the cost function
#include<config.h>
#include<utils/macros.h>
#include<fit/fit_wrapper/costfunction_policy.h>

#include<minimizer.h>

CPSFIT_START_NAMESPACE

#define INHERIT_MINIMIZER_POLICY_TYPEDEFS(FROM) \
  INHERIT_COSTFUNCTION_POLICY_TYPEDEFS(FROM);\
  INHERIT_TYPEDEF(FROM,minimizerType); \
  INHERIT_TYPEDEF(FROM,minimizerParamsType)


template<typename CostFuncPolicy>
class MarquardtLevenbergMinimizerPolicy: public CostFuncPolicy{
public:
  INHERIT_COSTFUNCTION_POLICY_TYPEDEFS(CostFuncPolicy);

  typedef MarquardtLevenbergMinimizer<costFunctionType> minimizerType;
  typedef typename minimizerType::AlgorithmParameterType minimizerParamsType;
};

template<typename CostFuncPolicy>
class GSLtrsMinimizerPolicy: public CostFuncPolicy{
public:
  INHERIT_COSTFUNCTION_POLICY_TYPEDEFS(CostFuncPolicy);

  typedef GSLtrsMinimizer<costFunctionType> minimizerType;
  typedef typename minimizerType::AlgorithmParameterType minimizerParamsType;
};


CPSFIT_END_NAMESPACE

#endif
