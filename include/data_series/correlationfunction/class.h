#ifndef _CPSFIT_CORRELATION_FUNCTION_CLASS_H_
#define _CPSFIT_CORRELATION_FUNCTION_CLASS_H_

#include<config.h>
#include<utils/macros.h>
#include<data_series/correlationfunction/types.h>
#include<data_series/data_series.h>

CPSFIT_START_NAMESPACE

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

CPSFIT_END_NAMESPACE
#endif
