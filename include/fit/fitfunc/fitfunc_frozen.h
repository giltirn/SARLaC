#ifndef _SARLAC_FITFUNC_FROZEN_H_
#define _SARLAC_FITFUNC_FROZEN_H_

//A "frozen fit" in which a fitfunc has some of its parameters fixed and unvaried during the minimization, can be implemented simply as a wrapper around a standard fit function
#include<config.h>
#include<utils/macros.h>
#include<containers/parameter_vector.h>

SARLAC_START_NAMESPACE

namespace _FrozenFitFunc_helper{
  template<typename T, bool default_constructible>
  struct _construct{};

  template<typename T>
  struct _construct<T, true>{
    inline static T construct(const std::unique_ptr<T> &freeze_vals){ return freeze_vals ? T(*freeze_vals) : T(); } //note the default constructor will not properly setup the output for variable-size types like NumericVector. For these types make sure you have provided a value for freeze_vals even if you are not going to freeze any parameters
  };
  template<typename T>
  struct _construct<T, false>{
    inline static T construct(const std::unique_ptr<T> &freeze_vals){
      if(!freeze_vals) error_exit(std::cout << "FrozenFitFunc::construct If the parameter type is not default constructible user must call 'freeze' with a value that can be copied (its contents do not matter)\n");
      return T(*freeze_vals);
    }
  };
};

//Implement frozen fits as wrapper around fit function with mapping. The object representing the 'superset', i.e. the parameters of the original fit function, must have an operator() returning a fixed type
template<typename FitFunc>
class FrozenFitFunc{
public:
  typedef typename std::remove_reference<decltype( ((typename FitFunc::ParameterType*)nullptr)->operator()(0) )>::type ParamElementType;
  typedef typename FitFunc::ParameterType ParameterSuperType;
  typedef typename FitFunc::ValueDerivativeType ValueDerivativeSuperType;

  typedef typename FitFunc::ValueType ValueType;
  typedef typename FitFunc::GeneralizedCoordinate GeneralizedCoordinate;
  typedef parameterVector<ParamElementType> ParameterType;
  typedef parameterVector<ParamElementType> ValueDerivativeType;
private:

  const FitFunc &fitfunc;
  std::vector<bool> param_freeze;
  std::unique_ptr<ParameterSuperType> freeze_vals;
  int n_frozen;

  //-1 for indices not in the subset
  std::vector<int> superset_subset_map;

public:
  FrozenFitFunc(const FitFunc &_fitfunc): fitfunc(_fitfunc), param_freeze(fitfunc.Nparams(),false), n_frozen(0), 
					  superset_subset_map(fitfunc.Nparams()){}

  ParameterSuperType mapParamsSubsetToSuperset(const ParameterType &params_subset) const{
    if(n_frozen != 0) assert(freeze_vals);
    ParameterSuperType superset = _FrozenFitFunc_helper::_construct<ParameterSuperType,std::is_default_constructible<ParameterSuperType>::value >::construct(freeze_vals);
    
    int subset_idx = 0;
    for(int i=0;i<fitfunc.Nparams();i++){
      if(!param_freeze[i])
	superset(i) = params_subset[subset_idx++];
      else
	superset(i) = (*freeze_vals)(i);
    }
    return superset;
  }
  ParameterType mapParamsSupersetToSubset(const ParameterSuperType &params_superset) const{
    ParameterType subset(Nparams());
    int subset_idx = 0;
    for(int i=0;i<fitfunc.Nparams();i++){
      if(!param_freeze[i]) subset[subset_idx++] = params_superset(i);
    }
    return subset;
  }

  //Returns -1 if parameter is not in the subset
  inline int getParamsSubsetIndex(const int superset_idx) const{
    return superset_subset_map[superset_idx];
  }

  //This should only be called once because subsequent calls overwrite the existing freeze information
  void freeze(const std::vector<int> &params, const ParameterSuperType &from){
    //Reset existing freeze info
    for(int i=0;i<param_freeze.size();i++) param_freeze[i] = false;  //reset
    
    freeze_vals.reset(new ParameterSuperType(from));
    for(int i=0;i<params.size();i++) param_freeze[params[i]] = true;
    n_frozen = params.size();
    
    int subset_idx = 0;
    for(int i=0;i<fitfunc.Nparams();i++)
      superset_subset_map[i] = param_freeze[i] ? -1 : subset_idx++;
  }
  
  inline ValueType value(const GeneralizedCoordinate &coord, const ParameterType &params_subset) const{
    return fitfunc.value(coord, mapParamsSubsetToSuperset(params_subset));
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &coord, const ParameterType &params_subset) const{
    ValueDerivativeType derivs_subset(Nparams());
    ValueDerivativeSuperType derivs_superset = fitfunc.parameterDerivatives(coord, mapParamsSubsetToSuperset(params_subset));
    int subset_idx = 0;
    for(int i=0;i<fitfunc.Nparams();i++)
      if(!param_freeze[i])
	derivs_subset[subset_idx++] = derivs_superset(i);
    return derivs_subset;
  }
  
  inline int Nparams() const{ return fitfunc.Nparams()-n_frozen; }
};

SARLAC_END_NAMESPACE
#endif
