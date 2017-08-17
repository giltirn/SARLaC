#ifndef _GENERAL_CONTAINER_H_
#define _GENERAL_CONTAINER_H_

//A way to store a generic objects in a common container
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

  inline bool is_null() const{ return !bool(v); }
  
  template<typename T>
  inline bool is() const{ assert(!is_null()); return v->is<T>(); }

  template<typename T>
  inline const T & value() const{ assert(!is_null()); return v->cast<T>()(); }

  template<typename T>
  inline T & value(){ assert(!is_null()); return v->cast<T>()(); }
};


#endif
