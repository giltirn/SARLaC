#ifndef _CPSFIT_TEMPLATE_WIZARDRY_DISTRIBUTION_H_
#define _CPSFIT_TEMPLATE_WIZARDRY_DISTRIBUTION_H_

//Metaprogramming constructs for obtaining information about distribution-like objects
#include<cstdlib>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>

CPSFIT_START_NAMESPACE

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


CPSFIT_END_NAMESPACE
#endif
