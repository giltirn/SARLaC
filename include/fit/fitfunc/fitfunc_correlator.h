#ifndef _SARLAC_FITFUNC_CORRELATOR_H_
#define _SARLAC_FITFUNC_CORRELATOR_H_

//Fit functions for lattice correlation functions: A*cosh(m*t), A*sinh(m*t), A*exp(-m*t)

#include<cmath>

#include<config.h>
#include<utils/macros.h>
#include<parser/parser.h>
#include<ET/generic_ET.h>
#include<containers/enumerated_struct.h>
#include<containers/parameter_vector.h>
#include<tensors/dual_number.h>

SARLAC_START_NAMESPACE

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

  ParameterType guess() const{ return ParameterType(1,0.5); }
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

  ParameterType guess() const{ return ParameterType(1,0.5); }
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

  ParameterType guess() const{ return ParameterType(1.0,0.5); }
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

  ParameterType guess() const{ return ParameterType(1,0.5,1,1); }
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

  ParameterType guess() const{ return ParameterType(1,0.5,1,1); }
};


//Fit N exponentials
//This also serves as an example of how to use dual numbers to automatically compute the derivative
class FitNStateExp{
public:
  int N;
  typedef double ValueType;
  typedef parameterVector<double> ParameterType; //format (A0,E0,A1,E1,...)
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord
  
  FitNStateExp(const int N): N(N){}

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate t, const ParameterType &p, const Boost &boost) const{
    if(p.size() != 2*N) error_exit(std::cout << "FitNStateExp::eval got parameter vector size " << p.size() << " expected " << 2*N << std::endl);
    T out(0);
    for(int j=0; j<2*N; j+=2){
      T A = boost(p[j],j);
      T E = boost(p[j+1],j+1);
      out = out + A*exp(-E*t);
    }
    return out;
  }
  inline ValueType value(const GeneralizedCoordinate t, const ParameterType &p) const{
    return eval<double>(t,p,
			[&](const double v, const int i){ return v; }
			);
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate t, const ParameterType &p) const{
    ValueDerivativeType d(2*N);
    for(int i=0;i<2*N;i++) d(i) = eval<dual>(t,p,
					     [&](const double v, const int j){ return dual(v, j==i ? 1.:0.); } //xp elem 1 for term we want to diff wrt
					     ).xp; //derivative component
    return d;
  }

  inline int Nparams() const{ return 2*N; }

};


SARLAC_END_NAMESPACE
#endif
