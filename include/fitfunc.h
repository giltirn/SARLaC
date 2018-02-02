#ifndef _FITFUNC_H
#define _FITFUNC_H

#include<numeric_tensors.h>
#include<template_wizardry.h>
#include<parser.h>

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type getCoord(const int i, const T coord){ return coord; }

template<typename T>
typename std::enable_if<is_std_vector<T>::value, typename T::value_type>::type getCoord(const int i, const T &coord){ return coord[i]; }

template<typename T, typename std::enable_if<!is_std_vector<T>::value,int>::type = 0 >
auto getCoord(const int i, const T &coord, ...)->decltype( coord[0] ){ return coord[i]; } //generic coordinate must have operator[]



//A multi-dimensional linear fit. 0-dimensional: Constant
template<typename _GeneralizedCoordinate, typename Numeric, int Dimension,
	 typename _ParameterType = NumericVector<Numeric>,  typename _ValueDerivativeType = NumericVector<Numeric> >  //p[0] + p[1]*x[0] + p[2]*x[1] + ...
class NumericLinearFit{
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


//Fit functions for 'standard fits'; exponential, cosh, sinh
struct StandardFitParams{
  double A;
  double m;
  StandardFitParams(){}
  StandardFitParams(const double _A, const double _m): A(_A), m(_m){}
  typedef StandardFitParams ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,StandardFitParams>::value, int>::type = 0>
  StandardFitParams(U&& expr){
    this->A = expr[0];
    this->m = expr[1];
  }
  inline double &operator()(const int i){ return i == 1 ? m : A; }
  inline const double &operator()(const int i) const{ return i == 1 ? m : A; }
  inline size_t size() const{ return 2;}

  inline void zero(){ A=m=0.; }
  inline std::string print() const{ std::ostringstream os; os << "A: " << A << " m: " << m; return os.str(); }
};
GENERATE_PARSER( StandardFitParams , (double,A)(double,m) );

template<>
struct getElem<StandardFitParams>{
  static inline auto elem(const StandardFitParams &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const StandardFitParams &v){ return 0; }
};


struct StandardFitParamDerivs{
  double dA;
  double dm;
  inline double &operator()(const int i){ return i == 1 ? dm : dA; }
  inline const double &operator()(const int i) const{ return i == 1 ? dm : dA; }
  inline size_t size() const{ return 2;}
  inline void zero(){ dA=dm=0.; }
  inline std::string print() const{ std::ostringstream os; os << "dA: " << dA << " dm: " << dm; return os.str(); }
  inline void resize(const int sz){ assert(sz == 2); }
};

struct StandardFitFuncBase{
  typedef double ValueType;
  typedef StandardFitParams ParameterType;
  typedef StandardFitParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  virtual ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const = 0;
  virtual ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const = 0;
  virtual int Nparams() const = 0;
};


class FitCosh: public StandardFitFuncBase{
  const double Lt;
public:
  typedef double ValueType;
  typedef StandardFitParams ParameterType;
  typedef StandardFitParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitCosh(const double _Lt): Lt(_Lt){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * ( exp(-p.m*t) + exp(-p.m*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = exp(-p.m*t) + exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*exp(-p.m*t) + -(Lt-t)*exp(-p.m*(Lt-t)) );
    return yderivs;
  }

  inline int Nparams() const{ return 2; }

  ParameterType guess(){ return ParameterType(1,0.5); }
};
class FitSinh: public StandardFitFuncBase{
  const double Lt;
public:
  typedef double ValueType;
  typedef StandardFitParams ParameterType;
  typedef StandardFitParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitSinh(const double _Lt): Lt(_Lt){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * ( exp(-p.m*t) - exp(-p.m*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = exp(-p.m*t) - exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*exp(-p.m*t) - -(Lt-t)*exp(-p.m*(Lt-t)) );
    return yderivs;
  }

  inline int Nparams() const{ return 2; }

  ParameterType guess(){ return ParameterType(1,0.5); }
};
class FitExp: public StandardFitFuncBase{
public:
  typedef double ValueType;
  typedef StandardFitParams ParameterType;
  typedef StandardFitParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * exp(-p.m*t);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = exp(-p.m*t);
    yderivs.dm = -p.A * t* exp(-p.m*t);
    return yderivs;
  }

  inline int Nparams() const{ return 2; }

  ParameterType guess(){ return ParameterType(1.0,0.5); }
};

class FitConstant{
public:
  typedef double ValueType;
  typedef NumericVector<double> ParameterType;
  typedef NumericVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p(0);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return ValueDerivativeType(1,1.0);
  }

  inline int Nparams() const{ return 1; }

  ParameterType guess(){ return ParameterType(1,1.0); }
};

#endif
