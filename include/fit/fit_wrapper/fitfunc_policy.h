#ifndef _CPSFIT_FIT_WRAPPER_FITFUNC_POLICY_H_
#define _CPSFIT_FIT_WRAPPER_FITFUNC_POLICY_H_

//The fit func policy contains and manipulates the fit function. It stores/points to the fit function and implements the extraction of the function's parameters for a given sample. In this way it can wrap invisibly around the use of a frozen fit function

#include<config.h>
#include<utils/macros.h>
#include<fit/fit_wrapper/base_typedefs.h>
#include<fit/fitfunc/fitfunc_frozen.h>
#include<fit/param_bounds.h>

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
  typedef typename DistributionType::template rebase<typename baseFitFunc::ParameterType> FitParameterDistribution; //going to intercept and transform to internal representation

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


enum class ParameterBound { Min, Max, Window };

//Bounded fits
template<typename InputFitTypes, typename _fitFunc>
class boundedFitFuncPolicy: public baseFitTypedefs<InputFitTypes>{
public:
  INHERIT_BASE_FIT_TYPEDEFS(baseFitTypedefs<InputFitTypes>);

  struct transform{
    int param;
    ParameterBound bound;
    double min;
    double max;
    
    transform(int param, ParameterBound bound, double min, double max): param(param), bound(bound), min(min), max(max){}

    template<typename Vtype>
    inline void mapBoundedToUnbounded(Vtype &v) const{
      switch(bound){
      case ParameterBound::Min:
	v(param) = MinBoundedMapping::mapBoundedToUnbounded(v(param),min);
	break;
      case ParameterBound::Max:
	v(param) = MaxBoundedMapping::mapBoundedToUnbounded(v(param),max);
	break;
      case ParameterBound::Window:
	v(param) = WindowBoundedMapping::mapBoundedToUnbounded(v(param),min,max);
	break;
      }
    }
    template<typename Vtype>
    inline void mapUnboundedToBounded(Vtype &v) const{
      switch(bound){
      case ParameterBound::Min:
	v(param) = MinBoundedMapping::mapUnboundedToBounded(v(param),min);
	break;
      case ParameterBound::Max:
	v(param) = MaxBoundedMapping::mapUnboundedToBounded(v(param),max);
	break;
      case ParameterBound::Window:
	v(param) = WindowBoundedMapping::mapUnboundedToBounded(v(param),min,max);
	break;
      }
    }

    template<typename derivType, typename Vtype>
    inline void jacobian(derivType &dv, const Vtype &v) const{ //v is unbounded variable
      switch(bound){
      case ParameterBound::Min:
	dv(param) = dv(param) * MinBoundedMapping::derivBoundedWrtUnbounded(v(param),min);
	break;
      case ParameterBound::Max:
	dv(param) = dv(param) * MaxBoundedMapping::derivBoundedWrtUnbounded(v(param),max);
	break;
      case ParameterBound::Window:
	dv(param) = dv(param) * WindowBoundedMapping::derivBoundedWrtUnbounded(v(param),min,max);
	break;
      }
    }

  };
  struct boundFFintercept{
    _fitFunc const* ff;
    std::vector<transform> const *trans;
    
    boundFFintercept(_fitFunc const* ff, std::vector<transform> const *trans): ff(ff), trans(trans){}

    typedef typename _fitFunc::ValueType ValueType;
    typedef typename _fitFunc::ParameterType ParameterType;
    typedef typename _fitFunc::ValueDerivativeType ValueDerivativeType;
    typedef typename _fitFunc::GeneralizedCoordinate GeneralizedCoordinate; 

    inline ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{ 
      ParameterType pb(p);
      for(auto it=trans->begin(); it != trans->end(); ++it) it->mapUnboundedToBounded(pb);
      return ff->value(t,pb); 
    }

    inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
      ParameterType pb(p);
      for(auto it=trans->begin(); it != trans->end(); ++it) it->mapUnboundedToBounded(pb);
      ValueDerivativeType vd = ff->parameterDerivatives(t,pb);
      for(auto it=trans->begin(); it != trans->end(); ++it) it->jacobian(vd, p);
      return vd;
    }

    int Nparams() const{ return ff->Nparams(); }
  };


  typedef _fitFunc baseFitFunc;
  typedef boundFFintercept fitFunc;
  typedef typename DistributionType::template rebase<typename fitFunc::ParameterType> FitParameterDistribution;

  class fitFuncPolicyState{
    std::unique_ptr<fitFunc> ff_b;
    std::vector<transform> const *trans;
  public:
    fitFuncPolicyState(baseFitFunc const* ff, std::vector<transform> const *trans): trans(trans){
      assert(ff != NULL);
      ff_b.reset(new fitFunc(ff, trans));
    }
    fitFuncPolicyState(fitFuncPolicyState &&r) = default;
    
    const fitFunc & getFitFunc() const{
      return *ff_b;
    }
    //Map input guess from bounded range to internal, unbounded range 
    typename fitFunc::ParameterType getParamsSample(const FitParameterDistribution &params, const int s){ 
      typename fitFunc::ParameterType out(iterate<FitParameterDistribution>::at(s,params));
      for(auto it=trans->begin(); it != trans->end(); ++it) it->mapBoundedToUnbounded(out);
      return out;
    }
    //Performed after mimimization complete, transform from unbounded to bounded range
    void setParamsSample(FitParameterDistribution &params, const typename fitFunc::ParameterType &p_s, const int s){ 
      auto &into = iterate<FitParameterDistribution>::at(s,params);
      into = p_s;
      for(auto it=trans->begin(); it != trans->end(); ++it) it->mapUnboundedToBounded(into);
    }    
  };
private:
  baseFitFunc const* ff;
  std::vector<transform> trans;
public:
  inline boundedFitFuncPolicy(): ff(NULL){}
  
  inline void importFitFunc(const baseFitFunc &fitfunc){ ff = &fitfunc; }

  inline void setBound(int param, const ParameterBound bound, double min, double max){ //min or max will be ignored as appropriate if using Max or Min bounds
    trans.push_back( transform(param, bound, min, max) );
  }

  inline fitFuncPolicyState generateFitFuncPolicyState(const int s){
    return fitFuncPolicyState(ff,&trans);
  }
};



CPSFIT_END_NAMESPACE
#endif
