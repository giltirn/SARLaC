#ifndef _FITFUNC_MAPPING_H
#define _FITFUNC_MAPPING_H

//A generic wrapper for a fitfunction allowing mapping between different parameter types. Can be used to implement frozen fits for example (although the actual method I use is more targeted and less complex). It may also be used for "global fits" where parameters of component fit functions are mapped to each other.

#include<type_traits>
#include<vector>
#include<memory>
#include<iostream>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<fit/fitfunc/fitfunc_mapping_elems.h>

SARLAC_START_NAMESPACE

template<typename FitFunc,
	 typename _ParameterSubsetType = typename FitFunc::ParameterType,
	 typename _ValueDerivativeSubsetType = typename FitFunc::ValueDerivativeType
	 >
class FitFuncGenericRemap{
public:  
  typedef _ParameterSubsetType ParameterSubsetType;
  typedef _ValueDerivativeSubsetType ValueDerivativeSubsetType;
  typedef typename FitFunc::ParameterType ParameterSupersetType;
  typedef typename FitFunc::ValueDerivativeType ValueDerivativeSupersetType;

  //Used by minimizer and whatnot as exposed interface
  typedef typename FitFunc::ValueType ValueType;
  typedef typename FitFunc::GeneralizedCoordinate GeneralizedCoordinate;
  typedef ParameterSubsetType ParameterType;
  typedef ValueDerivativeSubsetType ValueDerivativeType;

private:
  std::vector<ParameterMapBaseWrap<ParameterSupersetType, ParameterSubsetType> > param_map; //map from reduced -> full
  std::vector<ParameterMapBaseWrap<ValueDerivativeSubsetType, ValueDerivativeSupersetType> > deriv_map; //map from full -> reduced

  const FitFunc &fitfunc;
  
  //Defaults
  ParameterSupersetType param_superset_default;
  ValueDerivativeSubsetType deriv_subset_default;

  //Apply mappings
  inline ParameterSupersetType mapParamsSubsetToSuperset(const ParameterSubsetType &subset) const{
    ParameterSupersetType out(param_superset_default);
    for(int i=0;i<param_map.size();i++)
      param_map[i].copy(out, subset);   
    return out;
  }
  
  inline ValueDerivativeSubsetType mapDerivativesSupersetToSubset(const ValueDerivativeSupersetType &derivs_superset) const{
    ValueDerivativeSubsetType derivs_subset = deriv_subset_default;
    for(int i=0;i<deriv_map.size();i++)
      deriv_map[i].copy(derivs_subset, derivs_superset);
    return derivs_subset;
  }
  
public:
  FitFuncGenericRemap(const FitFunc &_fitfunc, const ParameterSupersetType _param_superset_default, const ValueDerivativeSubsetType &_deriv_subset_default):
    fitfunc(_fitfunc), param_superset_default(_param_superset_default), deriv_subset_default(_deriv_subset_default){}

  //Set mappings
  inline void addParamMapSubsetToSuperset(const ParameterMapBase<ParameterSupersetType, ParameterSubsetType> &param_map_subset_to_superset){
    param_map.push_back(ParameterMapBaseWrap<ParameterSupersetType, ParameterSubsetType>(param_map_subset_to_superset));
  }
  inline void addDerivMapSupersetToSubset(const ParameterMapBase<ValueDerivativeSubsetType, ValueDerivativeSupersetType> &deriv_map_superset_to_subset){
    deriv_map.push_back(ParameterMapBaseWrap<ValueDerivativeSubsetType, ValueDerivativeSupersetType>(deriv_map_superset_to_subset));
  }

  inline ValueType value(const GeneralizedCoordinate &coord, const ParameterType &params_subset) const{
    return fitfunc.value(coord, mapParamsSubsetToSuperset(params_subset));
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &coord, const ParameterType &params_subset) const{
    return mapDerivativesSupersetToSubset(fitfunc.parameterDerivatives(coord, mapParamsSubsetToSuperset(params_subset)));
  }
  
  int Nparams() const{ return deriv_subset_default.size(); } //number of parameters in subset
};

SARLAC_END_NAMESPACE

#endif
