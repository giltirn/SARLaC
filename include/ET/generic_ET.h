#ifndef _GENERIC_ET_H
#define _GENERIC_ET_H

#include<ET/generic_ET/binop_vv.h>
#include<ET/generic_ET/binop_sv.h>
#include<ET/generic_ET/binop_vs.h>
#include<ET/generic_ET/unary.h>
#include<ET/generic_ET/loop_eval.h>

CPSFIT_START_NAMESPACE

//An expression template engine for any 'vector' class (i.e. any class with elements that can be iterated over)
//For classes without an operator[] or size() methods, you need to create a specialization of getElem (cf ET/generic_ET/getelem.h)

//By default it creates vector-vector binary operators  + - * /
//For your class you can disable one or more of these by specializing disableGenericETbinOp

//It also defines vector-scalar binary operators + - * /   and scalar-vector + - *  (not /)

//Finally it defines unary operators - sqrt and exp

//Put this inside your class to enable the ET
//Tag is used to discriminate between classes of object; a binary op requires both ops have the same tag

//Example  ENABLE_GENERIC_ET(MyClass, MyClass<T>, MyClass<T>)

#define ENABLE_GENERIC_ET(CLASS_NAME, CLASS_FULL, TAG)			\
  typedef TAG ET_tag;							\
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,CLASS_FULL>::value, int>::type = 0> \
  CLASS_NAME(U&& expr): CLASS_NAME(expr.common_properties()){			\
    loop_eval<CLASS_FULL, U>(*this, std::forward<U>(expr));		\
  }\
  \
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,CLASS_FULL>::value, int>::type = 0> \
  CLASS_FULL & operator=(U && expr){ \
    this->resize(expr.common_properties());    \
    loop_eval<CLASS_FULL, U>(*this, std::forward<U>(expr));	\
    return *this; \
  }

CPSFIT_END_NAMESPACE

#endif
