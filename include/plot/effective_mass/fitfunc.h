#ifndef _SARLAC_EFFMASS_FITFUNC_H_
#define _SARLAC_EFFMASS_FITFUNC_H_

#include<config.h>
#include<utils/macros.h>
#include<containers/single_value_container.h>

SARLAC_START_NAMESPACE

//For fit functions with a single amplitude coefficient A*f(m,t) we can extract the effective mass by numerically inverting   C(t)/C(t+1) = f(m,t)/f(m,t+1)
template<typename FitFunc>
class Fit2ptEffectiveMass{
public:
  typedef typename FitFunc::ParameterType BaseParameterType;
  typedef typename FitFunc::ValueDerivativeType BaseDerivativeType;
  typedef double ValueType;
  typedef singleValueContainer<double> ParameterType;
  typedef singleValueContainer<double> ValueDerivativeType;
  typedef double GeneralizedCoordinate;
private:  
  FitFunc const* fitfunc;
  int params_mass_index; //index of the parameter that we are trying to determine
  BaseParameterType base;
public:

  Fit2ptEffectiveMass(const FitFunc &ff, const BaseParameterType _base, const int pmi): fitfunc(&ff),  params_mass_index(pmi), base(_base){
    //assert(base.size() == 2);
  }

  inline double value(const double t, const ParameterType &params) const{    
    BaseParameterType p(base); p(params_mass_index) = *params;
    return fitfunc->value(t,p)/fitfunc->value(t+1,p);
  }
  
  ValueDerivativeType parameterDerivatives(const double t, const ParameterType &params) const{
    ValueDerivativeType yderivs;
    BaseParameterType p(base); p(params_mass_index) = *params;
    
    double value_t = fitfunc->value(t,p);
    BaseDerivativeType subderivs_t = fitfunc->parameterDerivatives(t,p);

    double value_tp1 = fitfunc->value(t+1,p);
    BaseDerivativeType subderivs_tp1 = fitfunc->parameterDerivatives(t+1,p);
 
    *yderivs = subderivs_t(params_mass_index)/value_tp1 - value_t/value_tp1/value_tp1 * subderivs_tp1(params_mass_index);
    return yderivs;
  }

  inline int Nparams() const{ return 1; }
};


SARLAC_END_NAMESPACE
#endif
