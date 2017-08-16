#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H

//correlationFunction is a time series built on dataSeries but which has an expression-template engine for algebraic manipulations.
//User can modify how the ETE acts upon the underlying elements of the time series by changing the pair type
#include<generic_ET.h>
#include<data_series.h>

template<typename T, typename Tag>
struct tagged{
  T value;
  inline operator T() const{ return value; }
  inline tagged(const T &v): value(v){}
  inline tagged(){}
};

template<typename A, typename B>
struct CorrFuncTaggedPair: public std::pair<A,B>{
  using std::pair<A,B>::pair;
};


template<typename GeneralizedCoordinate, typename DistributionType, template<typename,typename> class PairType = CorrFuncTaggedPair>
class correlationFunction: public dataSeries<GeneralizedCoordinate, DistributionType, PairType>{
  typedef dataSeries<GeneralizedCoordinate, DistributionType, PairType> Parent;
public:
  typedef typename Parent::ElementType ElementType;
  typedef correlationFunction<GeneralizedCoordinate,DistributionType,PairType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,correlationFunction<GeneralizedCoordinate,DistributionType,PairType> >::value, int>::type = 0>
  correlationFunction(U&& expr) : Parent(expr.common_properties()){
   for(int i=0;i<this->size();i++)
     (*this)[i] = expr[i];
  }  
  explicit correlationFunction(const int n): Parent(n){}
  correlationFunction() = default;
  correlationFunction(correlationFunction &&r) = default;
  correlationFunction(const correlationFunction &r) = default;
  template<typename Initializer>
  inline correlationFunction(const int n, const Initializer &initializer): Parent(n,initializer){}

  correlationFunction<GeneralizedCoordinate,DistributionType,PairType> & operator=(const correlationFunction<GeneralizedCoordinate,DistributionType,PairType> &) = default;
  correlationFunction<GeneralizedCoordinate,DistributionType,PairType> & operator=(correlationFunction<GeneralizedCoordinate,DistributionType,PairType> &&) = default;
};

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
#endif
