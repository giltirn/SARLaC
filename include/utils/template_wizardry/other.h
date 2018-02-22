#ifndef _CPSFIT_TEMPLATE_OTHER_H_
#define _CPSFIT_TEMPLATE_OTHER_H_

#include<cstdlib>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry/types.h>

CPSFIT_START_NAMESPACE

//Check if a class T has an assignment operator to a type EqualsType
template<typename T, typename EqualsType, typename U = void>
struct hasEqualsMethod{
  enum{ value = 0 };
};
template<typename T, typename EqualsType>
struct hasEqualsMethod<T, EqualsType, typename Void<decltype( ((T*)(NULL))->operator=( *((EqualsType*)(NULL)) ) )>::type>{
  enum{ value = 1 };
};


CPSFIT_END_NAMESPACE
#endif
