#ifndef CPSFIT_PARAM_BOUNDS_H
#define CPSFIT_PARAM_BOUNDS_H

#include<config.h>
#include<utils/macros.h>

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

CPSFIT_END_NAMESPACE

#endif
