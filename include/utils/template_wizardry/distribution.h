#ifndef _SARLAC_TEMPLATE_WIZARDRY_DISTRIBUTION_H_
#define _SARLAC_TEMPLATE_WIZARDRY_DISTRIBUTION_H_

//Metaprogramming constructs for obtaining information about distribution-like objects
#include<cstdlib>

#include<config.h>
#include<utils/macros.h>
#include "types.h"
#include "type_classify.h"

SARLAC_START_NAMESPACE

//Check if a class T contains a typedef by name DataType
template<class T, class Fallback = void>
struct hasDataType{ enum{ value = 0 }; };

template<class T>
struct hasDataType<T, typename Void<typename T::DataType>::type >{ enum{ value = 1 }; };

//Check if a class T has a method sample(int)
template<typename T, typename U = void>
struct hasSampleMethod{
  enum{ value = 0 };
};
template<typename T>
struct hasSampleMethod<T, typename Void<decltype( ((T*)(NULL))->sample(0) )>::type>{
  enum{ value = 1 };
};

//Check if a class T has a method mean()
template<typename T, typename U = void>
struct hasMeanMethod{
  enum{ value = 0 };
};
template<typename T>
struct hasMeanMethod<T, typename Void<decltype( ((T*)(NULL))->mean() )>::type>{
  enum{ value = 1 };
};

//Check if a class T has a method zero()
template<typename T, typename U = void>
struct hasZeroMethod{
  enum{ value = 0 };
};
template<typename T>
struct hasZeroMethod<T, typename Void<decltype( ((T*)(NULL))->zero() )>::type>{
  enum{ value = 1 };
};

#define ENABLE_IF_HASSAMPLEMETHOD(Type)   typename std::enable_if< \
	   hasSampleMethod<Type>::value \
	   , int>::type = 0

template<typename T, ENABLE_IF_HASSAMPLEMETHOD(T)>
struct getSampleType{
  typedef typename std::decay< decltype( ( (const T*)NULL )->sample(0) ) >::type type;
};


//For compound distribution types, get the underlying POD base type
struct _get_base_type{
  struct _dist_mark;
  struct _other_mark;

  template<typename T>
  using classify = typename TypeIfElse< hasSampleMethod<T>::value, _dist_mark, _other_mark >::type;

  template<typename T,  typename mark>
  struct _get{};

  template<typename T>
  struct _get<T,_other_mark>{
    typedef T type;
  };

  template<typename T>
  struct _get<T,_dist_mark>{
    typedef typename _get<typename T::DataType, classify<typename T::DataType> >::type type;
  };
};
  
template<typename T>
struct getBaseType{
  typedef typename _get_base_type::_get<T, _get_base_type::classify<T> >::type type;
};




SARLAC_END_NAMESPACE
#endif
