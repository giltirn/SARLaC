#ifndef CPSFIT_PARAM_BOUNDS_H
#define CPSFIT_PARAM_BOUNDS_H

#include<config.h>
#include<utils/macros.h>
#include<parser/parser.h>

CPSFIT_START_NAMESPACE

//Fit parameter bounding can be achieved by transformations of bounded functions onto continous unbounded functions
//Specifically the input parameter is assumed to be bound and we map onto an unbounded range in which the minimizer performs its minimization

struct WindowBoundedMapping{
  static inline double mapUnboundedToBounded(const double pu, const double min, const double max){ //unbounded range -> bounded range
    return min + (sin(pu)+1)*(max-min)/2.;
  }
  static inline double mapBoundedToUnbounded(const double pb, const double min, const double max){ //bounded range -> unbounded range
    return asin(2.*(pb - min)/(max-min) - 1);
  }
  //derivative of bounded variable wrt unbounded variable
  //If 'u' is the unbounded variable and 'b' the bounded, then  df/du = df/db db/du
  static inline double derivBoundedWrtUnbounded(const double pu, const double min, const double max){  //'pu' is unbounded variable
    return cos(pu)*(max-min)/2.;
  }
};
struct MinBoundedMapping{
  static inline double mapUnboundedToBounded(const double pu, const double min){
    return min -1.+ sqrt( pu*pu + 1 );
  }
  static inline double mapBoundedToUnbounded(const double pb, const double min){
    return sqrt( pow( pb - min + 1., 2) - 1. );
  }
  static inline double derivBoundedWrtUnbounded(const double pu, const double min){
    return pu/sqrt( pu*pu + 1 );
  }
};
struct MaxBoundedMapping{
  static inline double mapUnboundedToBounded(const double pu, const double max){
    return max +1.- sqrt( pu*pu + 1 );
  }
  static inline double mapBoundedToUnbounded(const double pb, const double max){
    return sqrt( pow(pb - max - 1, 2) - 1. );
  }
  static inline double derivBoundedWrtUnbounded(const double pu, const double max){
    return -pu/sqrt( pu*pu + 1 );
  }
};

GENERATE_ENUM_AND_PARSER(ParameterBound, (Min)(Max)(Window) ); 

struct boundedParameterTransform{
  int param;
  ParameterBound bound;
  double min;
  double max;
    
  boundedParameterTransform(): param(0), bound(ParameterBound::Min), min(0.), max(0.){}
  boundedParameterTransform(int param, ParameterBound bound, double min, double max): param(param), bound(bound), min(min), max(max){}

  template<typename Vtype>
  inline void mapBoundedToUnbounded(Vtype &v) const{
    switch(bound){
    case ParameterBound::Min:
      v(param) = MinBoundedMapping::mapBoundedToUnbounded(v(param),min);
      break;
    case ParameterBound::Max:
      v(param) = MaxBoundedMapping::mapBoundedToUnbounded(v(param),max);
      break;
    case ParameterBound::Window:
      v(param) = WindowBoundedMapping::mapBoundedToUnbounded(v(param),min,max);
      break;
    }
  }
  template<typename Vtype>
  inline void mapUnboundedToBounded(Vtype &v) const{
    switch(bound){
    case ParameterBound::Min:
      v(param) = MinBoundedMapping::mapUnboundedToBounded(v(param),min);
      break;
    case ParameterBound::Max:
      v(param) = MaxBoundedMapping::mapUnboundedToBounded(v(param),max);
      break;
    case ParameterBound::Window:
      v(param) = WindowBoundedMapping::mapUnboundedToBounded(v(param),min,max);
      break;
    }
  }

  template<typename derivType, typename Vtype>
  inline void jacobian(derivType &dv, const Vtype &v) const{ //v is unbounded variable
    switch(bound){
    case ParameterBound::Min:
      dv(param) = dv(param) * MinBoundedMapping::derivBoundedWrtUnbounded(v(param),min);
      break;
    case ParameterBound::Max:
      dv(param) = dv(param) * MaxBoundedMapping::derivBoundedWrtUnbounded(v(param),max);
      break;
    case ParameterBound::Window:
      dv(param) = dv(param) * WindowBoundedMapping::derivBoundedWrtUnbounded(v(param),min,max);
      break;
    }
  }

};
  

GENERATE_PARSER(boundedParameterTransform, (int, param)(ParameterBound, bound)(double, min)(double, max) );

CPSFIT_END_NAMESPACE

#endif
