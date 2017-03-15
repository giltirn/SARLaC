#ifndef _FITFUNC_H
#define _FITFUNC_H

#include<numeric_tensors.h>
#include<template_wizardry.h>

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type getCoord(const int i, const T coord){ return coord; }

template<typename T>
typename std::enable_if<is_std_vector<T>::value, typename T::value_type>::type getCoord(const int i, const T &coord){ return coord[i]; }

template<typename T, typename std::enable_if<!is_std_vector<T>::value,int>::type = 0 >
auto getCoord(const int i, const T &coord, ...)->decltype( coord[0] ){ return coord[i]; } //generic coordinate must have operator[]



//A multi-dimensional linear fit. 0-dimensional: Constant
template<typename _GeneralizedCoordinate, typename Numeric, int Dimension>  //p[0] + p[1]*x[0] + p[2]*x[1] + ...
class NumericLinearFit{
public:
  typedef Numeric ValueType;
  typedef NumericVector<Numeric> ParameterType;
  typedef NumericVector<Numeric> ValueDerivativeType; //derivative wrt parameters
  typedef _GeneralizedCoordinate GeneralizedCoordinate;

  ValueType value(const GeneralizedCoordinate &coord, const ParameterType &params) const{
    ValueType out = params(0);
    for(int i=1;i<=Dimension;i++)
      out += params(i) * getCoord(i-1,coord);    
    return out;
  }
  void parameterDerivatives(ValueDerivativeType &yderivs, const GeneralizedCoordinate &coord, const ParameterType &params) const{
    yderivs.resize(Nparams());
    yderivs(0) = 1.;
    for(int i=1;i<=Dimension;i++){
      yderivs(i) = getCoord(i-1,coord);
    }
  }

  inline int Nparams() const{ return Dimension+1; }
};

#endif
