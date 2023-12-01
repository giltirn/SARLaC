#ifndef _SARLAC_TEMPLATE_OTHER_H_
#define _SARLAC_TEMPLATE_OTHER_H_

#include<cstdlib>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry/types.h>

SARLAC_START_NAMESPACE

//Check if a class T has an assignment operator to a type EqualsType
template<typename T, typename EqualsType, typename U = void>
struct hasEqualsMethod{
  enum{ value = 0 };
};
template<typename T, typename EqualsType>
struct hasEqualsMethod<T, EqualsType, typename Void<decltype( ((T*)(NULL))->operator=( *((EqualsType*)(NULL)) ) )>::type>{
  enum{ value = 1 };
};

//Added to a class "INTO_TYPE"'s body, for a data member "MEM" this will define a function copy_MEM that copies a member of the same name and type from a different class if it has such a member,
//otherwise nothing is copied 
#define DEF_GENERIC_COPY(INTO_TYPE,MEM)			      \
    template<typename B> static inline void copy_##MEM(...){} \
    template<typename B, decltype(INTO_TYPE::MEM) B::* bB = &B::MEM> static inline void copy_##MEM(INTO_TYPE &into, const B& from){ into.MEM= from.MEM; }

SARLAC_END_NAMESPACE
#endif
