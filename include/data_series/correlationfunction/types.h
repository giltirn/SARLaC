#ifndef _SARLAC_CORRELATION_FUNCTION_TYPES_H_
#define _SARLAC_CORRELATION_FUNCTION_TYPES_H_

#include<utility>

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

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

SARLAC_END_NAMESPACE
#endif
