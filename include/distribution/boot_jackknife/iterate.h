#ifndef _BOOT_JACKKNIFE_DIST_ITERATE_H_
#define _BOOT_JACKKNIFE_DIST_ITERATE_H_

#include<config.h>
#include<distribution/boot_jackknife/class.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE

template<typename T, template<typename> class V>
struct iterate<bootJackknifeDistribution<T,V> >{
  typedef T type;

  static inline int size(const bootJackknifeDistribution<T,V> &from){ 
    return from.origEnsJackknife().size() + from.size()*from.sample(0).size();
  }
  static inline const T& at(int i, const bootJackknifeDistribution<T,V> &from){
    int nn = from.origEnsJackknife().size();

    if(i < nn) 
      return from.origEnsJackknife().sample(i);
    i -= nn;

    nn = from.sample(0).size();
    
    return from.sample(i/nn).sample(i%nn);
  }
  static inline T & at(int i, bootJackknifeDistribution<T,V> &from){
    int nn = from.origEnsJackknife().size();

    if(i < nn) 
      return from.origEnsJackknife().sample(i);
    i -= nn;

    nn = from.sample(0).size();

    return from.sample(i/nn).sample(i%nn);
  }
  static inline std::vector<int> unmap(const int i, const bootJackknifeDistribution<T,V> &from){ 
    const int nn = from.sample(0).size();
    return std::vector<int>({i/nn, i%nn});
  }
};


template<typename T, template<typename> class V, int is_const>
class _distributionIterator<bootJackknifeDistribution<T,V>, is_const>{
  typedef typename add_const_if_int<bootJackknifeDistribution<T,V>,is_const>::type distributionType;
  typedef typename add_const_if_int<T,is_const>::type type;
  typedef iterate<bootJackknifeDistribution<T,V> > siter;
  int i;
  int sz;
  distributionType* d;
public:
  _distributionIterator(): i(0), sz(0){}
  _distributionIterator(distributionType &dist): i(0), sz(siter::size(dist)), d(&dist){}
  inline void operator++(){
    ++i;
  }
  inline bool end() const{ return i>=sz; }
  inline type& operator*() const{ return siter::at(i, *d); }
  inline void report() const{ std::cout << i << " (" << sz << ")" << std::endl; }
};



CPSFIT_END_NAMESPACE
#endif
