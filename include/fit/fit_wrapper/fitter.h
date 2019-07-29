#ifndef _CPSFIT_FIT_WRAPPER_FITTER_H_
#define _CPSFIT_FIT_WRAPPER_FITTER_H_

//Armed with an appropriate policy, the usual boilerplate fitting is implemented here

#include<config.h>
#include<utils/macros.h>
#include<fit/fit_wrapper/minimizer_policy.h>

CPSFIT_START_NAMESPACE

template<typename FitPolicies>
class fitter: public FitPolicies{
public:
  INHERIT_MINIMIZER_POLICY_TYPEDEFS(FitPolicies);
private:
  minimizerParamsType min_params;
public:
  fitter(){
    min_params.verbose = true;
  }
  fitter(const minimizerParamsType &_min_params): min_params(_min_params){}

  void setMinimizerParams(const minimizerParamsType &_min_params){ min_params = _min_params; }
  
  void fit(FitParameterDistribution &params,
	   DistributionType &chisq,
	   DistributionType &chisq_per_dof,
	   int &dof,
	   const CorrelationFunctionDistribution &data){

    const int ndata_fit = data.size();
    assert(ndata_fit > 0);
    const int nsample = data.value(0).size();
    assert(params.size() == nsample);
        
    chisq.resize(params.getInitializer());
    chisq_per_dof.resize(params.getInitializer());
    
    if(min_params.verbose) std::cout << "Starting fit with guess " << params << std::endl;

    int nosample = iterate<DistributionType>::size(data.value(0)); //some distributions have extra 'samples', eg jackknifeCdistribution has a sample representing the central value
    
    int nparams;
#pragma omp parallel for
    for(int s=0;s<nosample;s++){
      costFunctionState state = this->generateCostFunctionState(data,s);
      const costFunctionType &cost_func = state.getCostFunction();
      minimizerType minimizer(cost_func,min_params);

      typename fitFunc::ParameterType p_s = state.getParamsSample(params,s); //allows for arbitrary mapping from external to internal parameter type
      iterate<DistributionType>::at(s,chisq) = minimizer.fit(p_s);
      state.setParamsSample(params, p_s, s);
      
      iterate<DistributionType>::at(s,chisq_per_dof) = iterate<DistributionType>::at(s,chisq)/state.Ndof();
      assert(minimizer.hasConverged());
      if(s==0) dof = state.Ndof();
    }
  }

  inline void fit(FitParameterDistribution &params,
		  DistributionType &chisq,
		  DistributionType &chisq_per_dof,
		  const CorrelationFunctionDistribution &data){
    int dof;
    fit(params,chisq,chisq_per_dof,dof,data);
  }
};


CPSFIT_END_NAMESPACE
#endif
