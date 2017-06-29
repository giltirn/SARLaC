#ifndef PIPI_CORRELATION_FUNCTION_H
#define PIPI_CORRELATION_FUNCTION_H



template<typename T, typename Tag>
struct tagged{
  T value;
  inline operator T() const{ return value; }
  inline tagged(const T &v): value(v){}
  inline tagged(){}
};

template<typename A, typename B>
struct FitPiPiTaggedPair: public std::pair<A,B>{
  using std::pair<A,B>::pair;
};


template<typename DistributionType>
class correlationFunction: public dataSeries<double, DistributionType, FitPiPiTaggedPair>{
  typedef dataSeries<double, DistributionType, FitPiPiTaggedPair> Parent;
public:
  typedef typename Parent::ElementType ElementType;
  typedef correlationFunction<DistributionType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,correlationFunction<DistributionType> >::value, int>::type = 0>
  correlationFunction(U&& expr) : Parent(expr.common_properties()){
   for(int i=0;i<this->size();i++)
     (*this)[i] = expr[i];
  }
  explicit correlationFunction(const int n): Parent(n){}
  correlationFunction(correlationFunction &&r) = default;
  correlationFunction(const correlationFunction &r) = default;
  template<typename Initializer>
  inline correlationFunction(const int n, const Initializer &initializer): Parent(n,initializer){}

  correlationFunction<DistributionType> & operator=(const correlationFunction<DistributionType> &) = default;
  correlationFunction<DistributionType> & operator=(correlationFunction<DistributionType> &&) = default;
};

template<typename DistributionType>
struct getElem<correlationFunction<DistributionType> >{
  static inline auto elem(correlationFunction<DistributionType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline auto elem(const correlationFunction<DistributionType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline int common_properties(const correlationFunction<DistributionType> &v){
    return v.size();
  }
};

template<typename Dist>
using CFDpair = FitPiPiTaggedPair<double, Dist>;

template<typename Dist>
inline CFDpair<Dist> operator*(const int a, const CFDpair<Dist> &e){
  return CFDpair<Dist>(e.first, a*e.second);
}
template<typename Dist>
inline CFDpair<Dist> operator+(const CFDpair<Dist> &d, const CFDpair<Dist> &e){
  assert(e.first == d.first);
  return CFDpair<Dist>(e.first, d.second+e.second);
}
template<typename Dist>
inline CFDpair<Dist> operator-(const CFDpair<Dist> &d, const CFDpair<Dist> &e){
  assert(e.first == d.first);
  return CFDpair<Dist>(e.first, d.second-e.second);
}

#endif
