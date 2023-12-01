#ifndef _SARLAC_TEMPLATE_WIZARDRY_CONST_H_
#define _SARLAC_TEMPLATE_WIZARDRY_CONST_H_

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

//Add the const appelation to a type T if type U is const
template<typename T, typename U>
struct add_const_if{
  typedef T type;
};
template<typename T, typename U>
struct add_const_if<T, const U>{
  typedef const T type;
};

//Add const appelation to type T if add_const == 1
template<typename T, int add_const>
struct add_const_if_int{
  typedef T type;
};
template<typename T>
struct add_const_if_int<T,1>{
  typedef const T type;
};

SARLAC_END_NAMESPACE
#endif
