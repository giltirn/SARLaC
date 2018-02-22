#ifndef _CPSFIT_FIT_WRAPPER_POLICY_COMPOSE_H_
#define _CPSFIT_FIT_WRAPPER_POLICY_COMPOSE_H_

//Ultimately the policy for the wrapper comprises 3 layers: the base typedefs, the fitfunc policy and the cost func policy. To aid putting these together the user can obtain the type with the following compositor

#include<config.h>
#include<utils/macros.h>
#include<fit/fit_wrapper/base_typedefs.h>
#include<fit/fit_wrapper/fitfunc_policy.h>
#include<fit/fit_wrapper/costfunction_policy.h>

CPSFIT_START_NAMESPACE

template<typename FitFunc, template<typename,typename> class FitFuncPolicy, template<typename> class CostFunctionPolicy, typename InputFitTypes = standardInputFitTypes<typename FitFunc::GeneralizedCoordinate> >
struct composeFitPolicy{
  typedef CostFunctionPolicy<FitFuncPolicy<InputFitTypes,FitFunc> > type;
};

CPSFIT_END_NAMESPACE
#endif
