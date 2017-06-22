#ifndef PIPI_CORRELATION_FUNCTION_H
#define PIPI_CORRELATION_FUNCTION_H



template<typename T, typename Tag>
struct tagged{
  T value;
  inline operator T() const{ return value; }
  inline tagged(const T &v): value(v){}
  inline tagged(){}
};

template<typename DistributionType>
class correlationFunction: public dataSeries< tagged<double,correlationFunction<DistributionType> > , DistributionType>{
  typedef dataSeries< tagged<double,correlationFunction<DistributionType> > , DistributionType> Parent;
public:
  typedef typename Parent::ElementType ElementType;
  typedef correlationFunction<DistributionType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,correlationFunction<DistributionType> >::value, int>::type = 0>
  correlationFunction(U&& expr) : Parent(expr.common_properties()){
   for(int i=0;i<this->size();i++)
     (*this)[i] = expr[i];
  }

  correlationFunction(correlationFunction &&r) = default;
  template<typename Initializer>
  inline correlationFunction(const int n, const Initializer &initializer): Parent(n,initializer){}
};

template<typename DistributionType>
struct getElem<correlationFunction<DistributionType> >{
  static inline auto elem(correlationFunction<DistributionType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline auto elem(const correlationFunction<DistributionType> &v, const int i)->decltype(v[0]){ return v[i]; }
  static inline int common_properties(const correlationFunction<DistributionType> &v){
    return v.size();
  }
};

typedef std::pair<tagged<double, correlationFunction<distribution<double> > >, distribution<double> >  CFDpair;

inline CFDpair operator*(const int a, const CFDpair &e){
  return CFDpair(e.first, a*e.second);
}
inline CFDpair operator+(const CFDpair &d, const CFDpair &e){
  assert(e.first == d.first);
  return CFDpair(e.first, d.second+e.second);
}
inline CFDpair operator-(const CFDpair &d, const CFDpair &e){
  assert(e.first == d.first);
  return CFDpair(e.first, d.second-e.second);
}

#endif
