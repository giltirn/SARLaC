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
//DEF_ENUMERATED_STRUCT( ( StandardFitParams, double, (A)(m), (0.0)(0.0) ) ); 

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
    return p.A * ( ::exp(-p.m*t) + ::exp(-p.m*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = ::exp(-p.m*t) + ::exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*::exp(-p.m*t) + -(Lt-t)*::exp(-p.m*(Lt-t)) );
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
    return p.A * ( ::exp(-p.m*t) - ::exp(-p.m*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = ::exp(-p.m*t) - ::exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*::exp(-p.m*t) - -(Lt-t)*::exp(-p.m*(Lt-t)) );
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
    return p.A * ::exp(-p.m*t);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = ::exp(-p.m*t);
    yderivs.dm = -p.A * t* ::exp(-p.m*t);
    return yderivs;
  }

  inline int Nparams() const{ return 2; }

  ParameterType guess(){ return ParameterType(1.0,0.5); }
};

CPSFIT_END_NAMESPACE
#endif
