#ifndef _CPSFIT_CORRELATION_FUNCTION_ET_H_
#define _CPSFIT_CORRELATION_FUNCTION_ET_H_

#include<config.h>
#include<utils/macros.h>
#include<ET/generic_ET.h>
#include<data_series/correlationfunction/class.h>

CPSFIT_START_NAMESPACE

template<typename GeneralizedCoordinate, typename DistributionType,template<typename,typename> class PairType>
struct getElem<correlationFunction<GeneralizedCoordinate,DistributionType,PairType> >{
  static inline auto elem(correlationFunction<GeneralizedCoordinate,DistributionType,PairType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline auto elem(const correlationFunction<GeneralizedCoordinate,DistributionType,PairType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline int common_properties(const correlationFunction<GeneralizedCoordinate,DistributionType,PairType> &v){
    return v.size();
  }
};

template<typename Coord, typename Dist>
using CFDpair = CorrFuncTaggedPair<Coord, Dist>;

template<typename Coord, typename Dist>
inline CFDpair<Coord,Dist> operator*(const int a, const CFDpair<Coord,Dist> &e){
  return CFDpair<Coord,Dist>(e.first, a*e.second);
}
template<typename Coord, typename Dist>
inline CFDpair<Coord,Dist> operator+(const CFDpair<Coord,Dist> &d, const CFDpair<Coord,Dist> &e){
  assert(e.first == d.first);
  return CFDpair<Coord,Dist>(e.first, d.second+e.second);
}
template<typename Coord, typename Dist>
inline CFDpair<Coord,Dist> operator-(const CFDpair<Coord,Dist> &d, const CFDpair<Coord,Dist> &e){
  assert(e.first == d.first);
  return CFDpair<Coord,Dist>(e.first, d.second-e.second);
}
template<typename Coord, typename Dist>
inline CFDpair<Coord,Dist> operator/(const CFDpair<Coord,Dist> &d, const CFDpair<Coord,Dist> &e){
  assert(e.first == d.first);
  return CFDpair<Coord,Dist>(e.first, d.second/e.second);
}
template<typename Coord, typename Dist>
inline CFDpair<Coord,Dist> operator/(const CFDpair<Coord,Dist> &d, const double e){
  return CFDpair<Coord,Dist>(d.first, d.second/e);
}

CPSFIT_END_NAMESPACE
#endif
