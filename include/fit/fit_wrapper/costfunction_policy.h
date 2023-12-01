#ifndef _SARLAC_FIT_WRAPPER_COSTFUNC_POLICY_H_
#define _SARLAC_FIT_WRAPPER_COSTFUNC_POLICY_H_

//The cost function policy is responsible for handling the cost function used in the minimizer. It inherits methods from the fit function policy.
#include<fit/fit_wrapper/costfunction_policy/uncorrelated_fit_policy.h>
#include<fit/fit_wrapper/costfunction_policy/correlated_fit_policy.h>
#include<fit/fit_wrapper/costfunction_policy/correlated_cov_fit_policy.h>
#include<fit/fit_wrapper/costfunction_policy/frozen_correlated_fit_policy.h>

#define INHERIT_COSTFUNCTION_POLICY_TYPEDEFS(FROM) \
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FROM); \
  INHERIT_TYPEDEF(FROM,costFunctionType); \
  INHERIT_TYPEDEF(FROM,costFunctionState)

#endif
