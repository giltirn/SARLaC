#ifndef _FITFUNC_MAPPING_H
#define _FITFUNC_MAPPING_H



//A generic wrapper for a fitfunction allowing mapping between different parameter types. Can be used to implement frozen fits for example.
template<typename DataStructTo, typename DataStructFrom>
struct ParameterMapBase{
  virtual void copy(DataStructTo &to, const DataStructFrom &from) const = 0;
  virtual ParameterMapBase<DataStructTo,DataStructFrom> *clone() const = 0;
};

template<typename DataStructTo, typename DataStructFrom>
struct ParameterMapBaseWrap{
  ParameterMapBase<DataStructTo,DataStructFrom> *e;
public:
  inline void copy(DataStructTo &to, const DataStructFrom &from) const{ return e->copy(to,from); }
  ParameterMapBaseWrap(): e(NULL){}  
  ParameterMapBaseWrap(const ParameterMapBase<DataStructTo,DataStructFrom> &v): e(v.clone()){}
  ParameterMapBaseWrap(const ParameterMapBaseWrap<DataStructTo,DataStructFrom> &v): e(v.e->clone()){}
  
  ParameterMapBaseWrap &operator=(const ParameterMapBase<DataStructTo,DataStructFrom> &v){
    if(e!=NULL) delete e;
    e = v.clone();
  }
  
  ~ParameterMapBaseWrap(){ if(e!=NULL) delete e; }
};


template<typename T, typename DataStructTo, typename DataStructFrom>
struct ParameterMap: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toptr;
  T DataStructFrom::* fromptr;
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to.*toptr = from.*fromptr;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ParameterMap<T,DataStructTo, DataStructFrom>(*this); }
};
template<typename T, typename DataStructTo, typename DataStructFrom>
struct ParameterFreeze: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toptr;
  T value;
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to.*toptr = value;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ParameterFreeze<T,DataStructTo, DataStructFrom>(*this); }
};   
template<typename T, typename DataStructTo, typename DataStructFrom>
struct MemberArrayElementMap: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toarrayptr;
  T DataStructFrom::* fromarrayptr;
  int to_elem;
  int from_elem;
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    (to.*toarrayptr)[to_elem] = (from.*fromarrayptr)[from_elem];
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new MemberArrayElementMap<T,DataStructTo, DataStructFrom>(*this); }
};
template<typename T, typename DataStructTo, typename DataStructFrom>
struct MemberArrayElementFreeze: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toarrayptr;
  int to_elem;
  T value;
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    (to.*toarrayptr)[to_elem] = value;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new MemberArrayElementFreeze<T,DataStructTo, DataStructFrom>(*this); }
};
template<typename DataStructTo, typename DataStructFrom>
struct ArrayElementMap: public ParameterMapBase<DataStructTo, DataStructFrom>{
  int to_elem;
  int from_elem;

  ArrayElementMap(const int _to_elem, const int _from_elem): to_elem(_to_elem), from_elem(_from_elem){}
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to[to_elem] = from[from_elem];
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ArrayElementMap<DataStructTo, DataStructFrom>(*this); }
};
template<typename T, typename DataStructTo, typename DataStructFrom>
struct ArrayElementFreeze: public ParameterMapBase<DataStructTo, DataStructFrom>{
  int to_elem;
  T value;

  ArrayElementFreeze(const int _to_elem, const T &to): to_elem(_to_elem), value(to){}
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to[to_elem] = value;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ArrayElementFreeze<T,DataStructTo, DataStructFrom>(*this); }
};


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



//Implement frozen fits as wrapper around fit function with mapping
//This version is only for fit functions with parameters and derivatives stored in numerical vectors
//The above general mapping can also be used to do freezing but in a more cumbersome manner
template<typename FitFunc>
class FrozenFitFunc{
public:
  typedef typename FitFunc::ValueType ValueType;
  typedef typename FitFunc::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename FitFunc::ParameterType ParameterType;
  typedef typename FitFunc::ValueDerivativeType ValueDerivativeType;

private:
  typedef typename std::remove_reference<decltype( ((ParameterType*)nullptr)->operator[](0) )>::type ParamElementType;
  const FitFunc &fitfunc;
  std::vector<bool> param_freeze;
  std::vector<ParamElementType> freeze_vals;
  int n_frozen;
  
public:
  FrozenFitFunc(const FitFunc &_fitfunc): fitfunc(_fitfunc), param_freeze(fitfunc.Nparams(),false), freeze_vals(fitfunc.Nparams()), n_frozen(0){}

  ParameterType mapParamsSubsetToSuperset(const ParameterType &params_subset) const{
    ParameterType superset(fitfunc.Nparams());
    int subset_idx = 0;
    for(int i=0;i<fitfunc.Nparams();i++){
      if(param_freeze[i])
	superset[i] = freeze_vals[i];
      else
	superset[i] = params_subset[subset_idx++];
    }
    return superset;
  }
  ParameterType mapParamsSupersetToSubset(const ParameterType &params_superset) const{
    ParameterType subset(Nparams());
    int subset_idx = 0;
    for(int i=0;i<fitfunc.Nparams();i++){
      if(!param_freeze[i]) subset[subset_idx++] = params_superset[i];
    }
    return subset;
  }
  
  void freeze(const int param, const ParamElementType &value){
    param_freeze[param] = true;
    freeze_vals[param] = value;
    ++n_frozen;
  }
  
  inline ValueType value(const GeneralizedCoordinate &coord, const ParameterType &params_subset) const{
    return fitfunc.value(coord, mapParamsSubsetToSuperset(params_subset));
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &coord, const ParameterType &params_subset) const{
    ValueDerivativeType derivs_subset(Nparams());
    ValueDerivativeType derivs_superset = fitfunc.parameterDerivatives(coord, mapParamsSubsetToSuperset(params_subset));
    int subset_idx = 0;
    for(int i=0;i<fitfunc.Nparams();i++)
      if(!param_freeze[i])
	derivs_subset[subset_idx++] = derivs_superset[i];
    return derivs_subset;
  }
  
  int Nparams() const{ return fitfunc.Nparams()-n_frozen; }
};


#endif
