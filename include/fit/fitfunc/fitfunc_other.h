#ifndef _CPSFIT_FITFUNC_OTHER_H_
#define _CPSFIT_FITFUNC_OTHER_H_

#include<config.h>
#include<utils/macros.h>
#include<containers/parameter_vector.h>

CPSFIT_START_NAMESPACE

class FitConstant{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p(0);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return ValueDerivativeType(1,1.0);
  }

  inline int Nparams() const{ return 1; }

  ParameterType guess(){ return ParameterType(1,1.0); }
};

CPSFIT_END_NAMESPACE
#endif
