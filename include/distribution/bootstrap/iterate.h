#ifndef _BOOTSTRAP_DIST_ITERATE_H_
#define _BOOTSTRAP_DIST_ITERATE_H_

#include<config.h>
#include<distribution/bootstrap/class.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE


template<typename T, template<typename> class V, int is_const>
class _distributionIterator<bootstrapDistribution<T,V>,is_const>{
  typedef typename add_const_if_int<bootstrapDistribution<T,V>,is_const>::type distributionType;
  typedef typename add_const_if_int<T,is_const>::type type;
  int s;
  int sz;
  distributionType* d;
public:
  _distributionIterator(): s(0), sz(0), d(NULL){}
  _distributionIterator(distributionType &dist): s(0), sz(dist.size()+1), d(&dist){}
  inline void operator++(){ ++s; }
  inline bool end() const{ return s>=sz; }
  inline type& operator*() const{ return s==0 ? d->best() : d->sample(s-1); }
  inline void report() const{ std::cout << s << " (" << sz << ")" << std::endl; }
};

template<typename T, template<typename> class V>
struct iterate<bootstrapDistribution<T,V> >{
  static inline int size(const bootstrapDistribution<T,V> &from){ return from.size()+1; } 
  static inline const T& at(const int i, const bootstrapDistribution<T,V> &from){
    return i==0 ? from.best() : from.sample(i-1);
  }
  static inline T & at(const int i, bootstrapDistribution<T,V> &from){
    return i==0 ? from.best() : from.sample(i-1);
  }
  static inline std::vector<int> unmap(const int i, const bootstrapDistribution<T,V> &from){ 
    return std::vector<int>({i});
  }
};

CPSFIT_END_NAMESPACE
#endif
