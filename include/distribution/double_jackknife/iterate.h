#ifndef _DOUBLE_JACKKNIFE_DIST_ITERATE_H_
#define _DOUBLE_JACKKNIFE_DIST_ITERATE_H_

#include<config.h>
#include<distribution/double_jackknife/class.h>
#include<distribution/distribution_iterate.h>

SARLAC_START_NAMESPACE

template<typename T, template<typename> class V, int is_const>
class _distributionIterator<doubleJackknifeDistribution<T,V>, is_const>{
  typedef typename add_const_if_int<doubleJackknifeDistribution<T,V>,is_const>::type distributionType;
  typedef typename add_const_if_int<T,is_const>::type type;
  int s1;
  int s2;
  int sz1;
  int sz2;
  distributionType* d;
public:
  _distributionIterator(): s1(0), sz1(0), s2(0), sz2(0){}
  _distributionIterator(distributionType &dist): s1(0), sz1(dist.size()), s2(0), sz2(dist.size()-1), d(&dist){}
  inline void operator++(){
    ++s2; if(s2==sz2){ s2 = 0; ++s1; }
  }
  inline bool end() const{ return s1>=sz1; }
  inline type& operator*() const{ return d->sample(s1).sample(s2); }
  inline void report() const{ std::cout << s1 << " " << s2 << " (" << sz1 << " " << sz2 << ")" << std::endl; }
};

template<typename T, template<typename> class V>
struct iterate<doubleJackknifeDistribution<T,V> >{
  typedef T type;

  static inline int size(const doubleJackknifeDistribution<T,V> &from){ return from.size() * (from.size()-1); } //j + (from.size()-1)*i
  static inline const T& at(const int i, const doubleJackknifeDistribution<T,V> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline T & at(const int i, doubleJackknifeDistribution<T,V> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline std::vector<int> unmap(const int i, const doubleJackknifeDistribution<T,V> &from){ 
    const int nn = from.size()-1;
    return std::vector<int>({i/nn, i%nn});
  }
};

SARLAC_END_NAMESPACE
#endif
