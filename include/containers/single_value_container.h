#ifndef _CPSFIT_SINGLE_VALUE_CONTAINER_H_
#define _CPSFIT_SINGLE_VALUE_CONTAINER_H_

//In some cases we wish to interface with parts of the library that expect container-like objects with various canonical methods, but we wish to use a simple POD type. This class wraps a POD type with a container-like interface complete with the usual expression template engine

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<ET/generic_ET.h>

CPSFIT_START_NAMESPACE

template<typename T>
class singleValueContainer{
  T t;
public:
  ENABLE_GENERIC_ET(singleValueContainer, singleValueContainer<T>, singleValueContainer<T>);
  singleValueContainer(): t(0.){}
  explicit singleValueContainer(const T _t): t(_t){}
  inline T & operator()(const int i){ assert(i==0); return t; }
  inline const T &operator()(const int i) const{ assert(i==0); return t; }
  inline int size() const{ return 1; }
  inline std::string print() const{ return anyToStr<T>(t); }
  inline void resize(const int i){ if(i!=1) error_exit(std::cout << printType<singleValueContainer<T> >() << " resize called with value " << i << " != 1\n"); }
  inline void zero(){ t=0.; }
  inline T& operator*(){ return t; }
  inline const T& operator*() const{ return t; }
};
template<typename T>
struct getElem<singleValueContainer<T> >{
  static inline T& elem(singleValueContainer<T> &v, const int i){ return *v; }
  static inline const T& elem(const singleValueContainer<T> &v, const int i){ return *v; }    
  inline static int common_properties(const singleValueContainer<T> &v){ return 1; }
};

CPSFIT_END_NAMESPACE
#endif
