#ifndef _GENERIC_ET_GET_ELEM_H_
#define _GENERIC_ET_GET_ELEM_H_

#include<cstddef>

#include<utils/macros.h>

SARLAC_START_NAMESPACE

template<typename T>
struct getElem{
  static inline auto elem(const T &v, const int i)->decltype(v[i]){ return v[i]; }
  static inline auto elem(T &v, const int i)->decltype(v[i]){ return v[i]; }
  static inline size_t common_properties(const T &v){ return v.size(); }
};

SARLAC_END_NAMESPACE

#endif
