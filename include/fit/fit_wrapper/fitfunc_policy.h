#ifndef _CPSFIT_FIT_WRAPPER_FITFUNC_POLICY_H_
#define _CPSFIT_FIT_WRAPPER_FITFUNC_POLICY_H_

//The fit func policy contains and manipulates the fit function. It stores/points to the fit function and implements the extraction of the function's parameters for a given sample. In this way it can wrap invisibly around the use of a frozen fit function

#include<config.h>
#include<utils/macros.h>
#include<fit/fit_wrapper/base_typedefs.h>
#include<fit/fitfunc/fitfunc_frozen.h>

CPSFIT_START_NAMESPACE

//User-implemented policies should define the following types
#define INHERIT_FITFUNC_POLICY_TYPEDEFS(FROM)		       \
  INHERIT_BASE_FIT_TYPEDEFS(FROM);			       \
  INHERIT_TYPEDEF(FROM,fitFunc); \
  INHERIT_TYPEDEF(FROM,baseFitFunc); \
  INHERIT_TYPEDEF(FROM,FitParameterDistribution); \
  INHERIT_TYPEDEF(FROM,fitFuncPolicyState)

//Standard, unfrozen fits
template<typename InputFitTypes, typename _fitFunc>
class standardFitFuncPolicy: public baseFitTypedefs<InputFitTypes>{
public:
  INHERIT_BASE_FIT_TYPEDEFS(baseFitTypedefs<InputFitTypes>);
  typedef _fitFunc baseFitFunc;
  typedef _fitFunc fitFunc;
  typedef typename DistributionType::template rebase<typename fitFunc::ParameterType> FitParameterDistribution;

  class fitFuncPolicyState{
    fitFunc const* ff;
  public:
    fitFuncPolicyState(fitFunc const* _ff): ff(_ff){}
    const fitFunc & getFitFunc() const{
      assert(ff != NULL);
      return *ff;
    }
    typename fitFunc::ParameterType getParamsSample(const FitParameterDistribution &params, const int s){ return iterate<FitParameterDistribution>::at(s,params); }
    void setParamsSample(FitParameterDistribution &params, const typename fitFunc::ParameterType &p_s, const int s){ iterate<FitParameterDistribution>::at(s,params) = p_s; }
  };
private:
  fitFunc const* ff;
public:
  standardFitFuncPolicy(): ff(NULL){}
  
  void importFitFunc(const fitFunc &fitfunc){ ff = &fitfunc; }
  
  inline fitFuncPolicyState generateFitFuncPolicyState(const int s){
    return fitFuncPolicyState(ff);
  }
};


//Frozen fits
template<typename InputFitTypes, typename _fitFunc>
class frozenFitFuncPolicy: public baseFitTypedefs<InputFitTypes>{
public:
  INHERIT_BASE_FIT_TYPEDEFS(baseFitTypedefs<InputFitTypes>);

  typedef _fitFunc baseFitFunc;
  typedef FrozenFitFunc<_fitFunc> fitFunc;
  typedef jackknifeDistribution<typename baseFitFunc::ParameterType> FitParameterDistribution; //going to intercept and transform to internal representation

  class fitFuncPolicyState{
    std::unique_ptr<fitFunc> ff_frz;
  public:
    fitFuncPolicyState(baseFitFunc const* ff, const std::vector<int> &freeze_params, const std::unique_ptr<FitParameterDistribution> &freeze_vals, const int s){
      assert(ff != NULL);
      ff_frz.reset(new fitFunc(*ff));
      if(freeze_params.size() > 0){
	assert(freeze_vals);
	ff_frz->freeze(freeze_params, iterate<FitParameterDistribution>::at(s,*freeze_vals));
      }
    }
    fitFuncPolicyState(fitFuncPolicyState &&r): ff_frz(std::move(r.ff_frz)){}
    
    const fitFunc & getFitFunc() const{
      return *ff_frz;
    }
    typename fitFunc::ParameterType getParamsSample(const FitParameterDistribution &params, const int s){ return ff_frz->mapParamsSupersetToSubset(iterate<FitParameterDistribution>::at(s,params)); }
    void setParamsSample(FitParameterDistribution &params, const typename fitFunc::ParameterType &p_s, const int s){ iterate<FitParameterDistribution>::at(s,params) = ff_frz->mapParamsSubsetToSuperset(p_s); }    
  };
private:
  baseFitFunc const* ff;
  std::unique_ptr<FitParameterDistribution> freeze_vals;
  std::vector<int> freeze_params;
public:
  frozenFitFuncPolicy(): ff(NULL){}
  
  void importFitFunc(const baseFitFunc &fitfunc){ ff = &fitfunc; }

  void freeze(const std::vector<int> &_freeze_params, const FitParameterDistribution &_freeze_vals){
    freeze_vals.reset(new FitParameterDistribution(_freeze_vals));   freeze_params = _freeze_params;
  }
  
  inline fitFuncPolicyState generateFitFuncPolicyState(const int s){
    return fitFuncPolicyState(ff,freeze_params,freeze_vals,s);
  }
};


CPSFIT_END_NAMESPACE
#endif
