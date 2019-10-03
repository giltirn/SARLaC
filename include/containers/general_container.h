#ifndef _GENERAL_CONTAINER_H_
#define _GENERAL_CONTAINER_H_

#include<memory>
#include<cassert>

#include<containers/general_container/general_container_i.h>

CPSFIT_START_NAMESPACE

//A polymorphic container type
class generalContainer{
  std::unique_ptr<generalContainer_i> v;
public:
  generalContainer(){}
  template<typename T>
  generalContainer(const T &_v): v(new generalContainerEntry_i<T>(_v)){}
  generalContainer(const generalContainer &r): v(r.v ? r.v->clone() : NULL){}
  generalContainer(generalContainer &&r): v(r.v ? std::move(r.v) : NULL){}

  generalContainer &operator=(const generalContainer &r){ v.reset(r.v ? r.v->clone() : NULL); return *this; }
  generalContainer &operator=(generalContainer &&r){ v.reset(r.v ? r.v.release() : NULL); return *this; }
  template<typename T>
  generalContainer &operator=(const T &_v){ v.reset(new generalContainerEntry_i<T>(_v)); return *this; }
  template<typename T>
  generalContainer &operator=(T &&_v){ v.reset(new generalContainerEntry_i<T>(std::move(_v))); return *this; }

  inline bool is_null() const{ return !bool(v); }
  
  template<typename T>
  inline bool is() const{ assert(!is_null()); return v->is<T>(); }

  template<typename T>
  inline const T & value() const{ assert(!is_null()); return v->cast<T>()(); }

  template<typename T>
  inline T & value(){ assert(!is_null()); return v->cast<T>()(); }
};


CPSFIT_END_NAMESPACE

#endif
