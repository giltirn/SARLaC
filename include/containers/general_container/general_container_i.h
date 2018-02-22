#ifndef _GENERAL_CONTAINER_I_H__
#define _GENERAL_CONTAINER_I_H__

#include<utility>
#include<cstddef>

#include<utils/macros.h>

CPSFIT_START_NAMESPACE

template<typename T>
class generalContainerEntry_i;

class generalContainer_i{
public:
  template<typename T>
  inline bool is() const{
    generalContainerEntry_i<T> const* p = dynamic_cast<generalContainerEntry_i<T> const*>(this);
    return p != NULL;
  }
  template<typename T>
  const generalContainerEntry_i<T> & cast() const{ return *dynamic_cast<generalContainerEntry_i<T> const*>(this); }
  
  template<typename T>
  generalContainerEntry_i<T> & cast(){ return *dynamic_cast<generalContainerEntry_i<T>*>(this); }

  virtual generalContainer_i* clone() const = 0;  
  virtual ~generalContainer_i(){}
};

template<typename T>
class generalContainerEntry_i: public generalContainer_i{
  T v;
public:
  generalContainerEntry_i() = default;
  generalContainerEntry_i(const T &_v): v(_v){}
  generalContainerEntry_i(T &&_v): v(std::move(_v)){}
    
  const T & operator() ()const{ return v; }
  T & operator() (){ return v; }

  generalContainer_i* clone() const{ return new generalContainerEntry_i(*this); }
};

CPSFIT_END_NAMESPACE

#endif
