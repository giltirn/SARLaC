#ifndef _CPSFIT_FITFUNC_CORRELATOR_H_
#define _CPSFIT_FITFUNC_CORRELATOR_H_

//Fit functions for lattice correlation functions: A*cosh(m*t), A*sinh(m*t), A*exp(-m*t)

#include<cmath>

#include<config.h>
#include<utils/macros.h>
#include<parser/parser.h>
#include<ET/generic_ET.h>
#include<containers/enumerated_struct.h>

CPSFIT_START_NAMESPACE

//Fit functions for 'standard fits'; exponential, cosh, sinh
DEF_ENUMERATED_STRUCT( ( StandardFitParams, double, (A)(m), (0.0)(0.0) ) ); 

struct StandardFitFuncBase{
  typedef double ValueType;
  typedef StandardFitParams ParameterType;
  typedef StandardFitParams ValueDerivativeType; //derivative wrt parameters
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
  typedef StandardFitParams ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitCosh(const double _Lt): Lt(_Lt){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * ( ::exp(-p.m*t) + ::exp(-p.m*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.A = ::exp(-p.m*t) + ::exp(-p.m*(Lt-t));
    yderivs.m = p.A * ( -t*::exp(-p.m*t) + -(Lt-t)*::exp(-p.m*(Lt-t)) );
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
  typedef StandardFitParams ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitSinh(const double _Lt): Lt(_Lt){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * ( ::exp(-p.m*t) - ::exp(-p.m*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.A = ::exp(-p.m*t) - ::exp(-p.m*(Lt-t));
    yderivs.m = p.A * ( -t*::exp(-p.m*t) - -(Lt-t)*::exp(-p.m*(Lt-t)) );
    return yderivs;
  }

  inline int Nparams() const{ return 2; }

  ParameterType guess(){ return ParameterType(1,0.5); }
};
class FitExp: public StandardFitFuncBase{
public:
  typedef double ValueType;
  typedef StandardFitParams ParameterType;
  typedef StandardFitParams ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * ::exp(-p.m*t);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.A = ::exp(-p.m*t);
    yderivs.m = -p.A * t* ::exp(-p.m*t);
    return yderivs;
  }

  inline int Nparams() const{ return 2; }

  ParameterType guess(){ return ParameterType(1.0,0.5); }
};


//Two-state fits

DEF_ENUMERATED_STRUCT( ( TwoStateFitParams, double, (A0)(E0)(A1)(E1), (0.0)(0.0)(0.0)(0.0) ) ); 

class FitTwoStateCosh{
  const double Lt;
public:
  typedef double ValueType;
  typedef TwoStateFitParams ParameterType;
  typedef TwoStateFitParams ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitTwoStateCosh(const double _Lt): Lt(_Lt){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A0 * ( ::exp(-p.E0*t) + ::exp(-p.E0*(double(Lt)-t)) ) + p.A1 * ( ::exp(-p.E1*t) + ::exp(-p.E1*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.A0 = ::exp(-p.E0*t) + ::exp(-p.E0*(Lt-t));
    yderivs.E0 = p.A0 * ( -t*::exp(-p.E0*t) + -(Lt-t)*::exp(-p.E0*(Lt-t)) );
    yderivs.A1 = ::exp(-p.E1*t) + ::exp(-p.E1*(Lt-t));
    yderivs.E1 = p.A1 * ( -t*::exp(-p.E1*t) + -(Lt-t)*::exp(-p.E1*(Lt-t)) );
    return yderivs;
  }

  inline int Nparams() const{ return 4; }

  ParameterType guess(){ return ParameterType(1,0.5,1,1); }
};


class FitTwoStateExp{
public:
  typedef double ValueType;
  typedef TwoStateFitParams ParameterType;
  typedef TwoStateFitParams ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A0 * ::exp(-p.E0*t) + p.A1 * ::exp(-p.E1*t);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.A0 = ::exp(-p.E0*t);
    yderivs.E0 = p.A0 * ( -t*::exp(-p.E0*t) );
    yderivs.A1 = ::exp(-p.E1*t);
    yderivs.E1 = p.A1 * ( -t*::exp(-p.E1*t) );
    return yderivs;
  }

  inline int Nparams() const{ return 4; }

  ParameterType guess(){ return ParameterType(1,0.5,1,1); }
};


CPSFIT_END_NAMESPACE
#endif
