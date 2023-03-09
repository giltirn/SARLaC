#ifndef _CPSFIT_SIMPLE_FIT_WRAPPER_FITFUNC_WRAPPER_H
#define _CPSFIT_SIMPLE_FIT_WRAPPER_FITFUNC_WRAPPER_H

#include<config.h>
#include<utils/macros.h>

#include<containers/general_container.h>
#include<containers/parameter_vector.h>

CPSFIT_START_NAMESPACE

//Fit functions used by the simple fit wrapper should all inherit from this class. Below I provide a derived class that wraps an existing fit function so the user does not need to modify the fit function class itself
struct genericFitFuncBase{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef generalContainer GeneralizedCoordinate;

  virtual ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const = 0;
  virtual ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const = 0;
  virtual int Nparams() const = 0;
  //virtual ParameterType guess() const = 0;
  virtual ~genericFitFuncBase(){}

  //Some useful functions

  //Copy between two vector-like objects
  template<typename T, typename U>
  inline static void copyOver(T &out, const U &in){
    for(int i=0;i<in.size();i++) out(i) = in(i);
  }
  //Translate the coordinate into the underlying type T
  template<typename T>
  static inline T getCoord(const GeneralizedCoordinate &x){ assert(x.is<T>()); return x.value<T>(); }

  template<typename T>
  static inline GeneralizedCoordinate getWrappedCoord(const T &x){ return generalContainer(x); } 

  template<typename T>
  static inline ParameterType getWrappedParams(const T &x){ 
    ParameterType p(x.size()); copyOver(p,x); return p;
  } 

};  

#define INHERIT_GENERIC_FITFUNC_BASE_TYPEDEFS		\
  typedef genericFitFuncBase::ValueType ValueType;	\
  typedef genericFitFuncBase::ParameterType ParameterType;		\
  typedef genericFitFuncBase::ValueDerivativeType ValueDerivativeType;	\
  typedef genericFitFuncBase::GeneralizedCoordinate GeneralizedCoordinate


//This wrapper can be used to enwrap any function that confirms to the general template for use in the simplified fitter
//Note: for parameter classes which do not have a constructor that takes only the number of parameters
//      you should use the second wrapper below
template<typename FitFunc>
struct simpleFitFuncWrapper: public genericFitFuncBase{
  INHERIT_GENERIC_FITFUNC_BASE_TYPEDEFS;
  FitFunc fitfunc;

  simpleFitFuncWrapper(const FitFunc &fitfunc): fitfunc(fitfunc){}
  
  ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    typename FitFunc::ParameterType pbase(p.size()); copyOver(pbase, p);
    return fitfunc.value(getCoord<typename FitFunc::GeneralizedCoordinate>(x), pbase);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    typename FitFunc::ParameterType pbase(p.size()); copyOver(pbase, p);
    auto derivs = fitfunc.parameterDerivatives(getCoord<typename FitFunc::GeneralizedCoordinate>(x),pbase);
    ValueDerivativeType out(fitfunc.Nparams()); copyOver(out, derivs);
    return out;
  } 

  int Nparams() const{ return fitfunc.Nparams(); }
};

//This wrapper works for any parameter type that has a copy constructor
template<typename FitFunc>
struct genericFitFuncWrapper: public genericFitFuncBase{
  INHERIT_GENERIC_FITFUNC_BASE_TYPEDEFS;
  typedef typename FitFunc::ParameterType BaseParameterType;
  FitFunc fitfunc;
  BaseParameterType psetup;

  //The 'psetup' argument should be an instance of the object containing the wrapped fit function's parameters
  //The contents of psetup don't matter, it simply needs to be sized correctly
  genericFitFuncWrapper(const FitFunc &fitfunc, const BaseParameterType &psetup): fitfunc(fitfunc), psetup(psetup){}
  
  ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    typename FitFunc::ParameterType pbase(psetup); copyOver(pbase, p);
    return fitfunc.value(getCoord<typename FitFunc::GeneralizedCoordinate>(x), pbase);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    typename FitFunc::ParameterType pbase(psetup); copyOver(pbase, p);
    auto derivs = fitfunc.parameterDerivatives(getCoord<typename FitFunc::GeneralizedCoordinate>(x),pbase);
    ValueDerivativeType out(fitfunc.Nparams()); copyOver(out, derivs);
    return out;
  } 

  int Nparams() const{ return fitfunc.Nparams(); }
};


CPSFIT_END_NAMESPACE
#endif
