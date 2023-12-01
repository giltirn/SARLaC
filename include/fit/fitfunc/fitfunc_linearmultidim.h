#ifndef _SARLAC_FITFUNC_LINEAR_MULTIDIM_H_
#define _SARLAC_FITFUNC_LINEAR_MULTIDIM_H_

//A linear fit in N dimensions: p[0] + p[1]*x[0] + p[2]*x[1] + ...

#include<config.h>
#include<utils/macros.h>
#include<containers/parameter_vector.h>

SARLAC_START_NAMESPACE

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type getCoord(const int i, const T coord){ return coord; }

template<typename T>
typename std::enable_if<is_std_vector<T>::value, typename T::value_type>::type getCoord(const int i, const T &coord){ return coord[i]; }

template<typename T, typename std::enable_if<!is_std_vector<T>::value,int>::type = 0 >
auto getCoord(const int i, const T &coord, ...)->decltype( coord[0] ){ return coord[i]; } //generic coordinate must have operator[]


//Dimension 0 just includes the constant term
template<typename _GeneralizedCoordinate, typename Numeric, int Dimension,
	 typename _ParameterType = parameterVector<Numeric>,  typename _ValueDerivativeType = parameterVector<Numeric> >  //p[0] + p[1]*x[0] + p[2]*x[1] + ...
class FitFuncLinearMultiDim{
public:
  typedef Numeric ValueType;
  typedef _ParameterType ParameterType;
  typedef _ValueDerivativeType ValueDerivativeType; //derivative wrt parameters
  typedef _GeneralizedCoordinate GeneralizedCoordinate;

  ValueType value(const GeneralizedCoordinate &coord, const ParameterType &params) const{
    ValueType out = params(0);
    for(int i=1;i<=Dimension;i++)
      out += params(i) * getCoord(i-1,coord);    
    return out;
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &coord, const ParameterType &params) const{
    ValueDerivativeType yderivs(Nparams());
    yderivs(0) = 1.;
    for(int i=1;i<=Dimension;i++){
      yderivs(i) = getCoord(i-1,coord);
    }
    return yderivs;
  }

  inline int Nparams() const{ return Dimension+1; }
};

SARLAC_END_NAMESPACE
#endif
