#ifndef _BOOT_JACKKNIFE_DIST_ITERATE_H_
#define _BOOT_JACKKNIFE_DIST_ITERATE_H_

#include<config.h>
#include<distribution/boot_jackknife/class.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE

template<typename T, template<typename> class V, int is_const>
class _distributionIterator<bootJackknifeDistribution<T,V>, is_const>{
  typedef typename add_const_if_int<bootJackknifeDistribution<T,V>,is_const>::type distributionType;
  typedef typename add_const_if_int<T,is_const>::type type;
  int s1;
  int s2;
  int sz1;
  int sz2;
  distributionType* d;
public:
  _distributionIterator(): s1(0), sz1(0), s2(0), sz2(0){}
  _distributionIterator(distributionType &dist): s1(0), sz1(dist.size()), s2(0), sz2(dist.size()>0 ? dist.sample(0).size() : 0), d(&dist){}
  inline void operator++(){
    ++s2; if(s2==sz2){ s2 = 0; ++s1; }
  }
  inline bool end() const{ return s1>=sz1; }
  inline type& operator*() const{ return d->sample(s1).sample(s2); }
  inline void report() const{ std::cout << s1 << " " << s2 << " (" << sz1 << " " << sz2 << ")" << std::endl; }
};

template<typename T, template<typename> class V>
struct iterate<bootJackknifeDistribution<T,V> >{
  typedef T type;

  static inline int size(const bootJackknifeDistribution<T,V> &from){ 
    return (1 + from.size()) * from.origEnsJackknife().size();
  }
  static inline const T& at(int i, const bootJackknifeDistribution<T,V> &from){
    const int nn = from.origEnsJackknife().size();

    if(i < nn) 
      return from.origEnsJackknife().sample(i);
    i -= nn;

    return from.sample(i/nn).sample(i%nn);
  }
  static inline T & at(int i, bootJackknifeDistribution<T,V> &from){
    const int nn = from.origEnsJackknife().size();

    if(i < nn) 
      return from.origEnsJackknife().sample(i);
    i -= nn;

    return from.sample(i/nn).sample(i%nn);
  }
  static inline std::vector<int> unmap(const int i, const bootJackknifeDistribution<T,V> &from){ 
    const int nn = from.sample(0).size();
    return std::vector<int>({i/nn, i%nn});
  }
};

CPSFIT_END_NAMESPACE
#endif
