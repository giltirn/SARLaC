#ifndef _DISTRIBUTION_ITERATE_H
#define _DISTRIBUTION_ITERATE_H

template<typename T, int add_const>
struct add_const_if_int{
  typedef T type;
};
template<typename T>
struct add_const_if_int<T,1>{
  typedef const T type;
};

//This type of iterator is bound to a distribution instance
template<typename T_unconst, int is_const>
class _distributionIterator{
  typedef typename add_const_if_int<T_unconst,is_const>::type distributionType;
  typedef typename add_const_if_int<typename T_unconst::DataType,is_const>::type type;
  int s;
  int sz;
  distributionType* d;
public:
  _distributionIterator(): s(0), sz(0), d(NULL){}
  _distributionIterator(distributionType &dist): s(0), sz(dist.size()), d(&dist){}
  inline void operator++(){ ++s; }
  inline bool end() const{ return s>=sz; }
  inline type& operator*() const{ return d->sample(s); }
  inline void report() const{ std::cout << s << " (" << sz << ")" << std::endl; }
};

template<typename T, int is_const>
class _distributionIterator<jackknifeCdistribution<T>,is_const>{
  typedef typename add_const_if_int<jackknifeCdistribution<T>,is_const>::type distributionType;
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

template<typename T, int is_const>
class _distributionIterator<doubleJackknifeDistribution<T>, is_const>{
  typedef typename add_const_if_int<doubleJackknifeDistribution<T>,is_const>::type distributionType;
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

template<typename T> using distributionIterator = _distributionIterator<typename std::remove_const<T>::type, std::is_const<T>::value>;





//These structures allow you to generically iterate over the elements of many distributions (of the same size)
template<typename distributionType>
struct iterate;

template<typename T>
struct iterate<doubleJackknifeDistribution<T> >{
  static inline int size(const doubleJackknifeDistribution<T> &from){ return from.size() * (from.size()-1); } //j + (from.size()-1)*i
  static inline const T& at(const int i, const doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline T & at(const int i, doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline std::vector<int> unmap(const int i, const doubleJackknifeDistribution<T> &from){ 
    const int nn = from.size()-1;
    return std::vector<int>({i/nn, i%nn});
  }
};
template<typename T>
struct iterate<rawDataDistribution<T> >{
  static inline int size(const rawDataDistribution<T> &from){ return from.size(); } 
  static inline const T& at(const int i, const rawDataDistribution<T> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, rawDataDistribution<T> &from){
    return from.sample(i);
  }
  static inline std::vector<int> unmap(const int i, const rawDataDistribution<T> &from){ 
    return std::vector<int>({i});
  }
};
template<typename T>
struct iterate<jackknifeDistribution<T> >{
  static inline int size(const jackknifeDistribution<T> &from){ return from.size(); } 
  static inline const T& at(const int i, const jackknifeDistribution<T> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, jackknifeDistribution<T> &from){
    return from.sample(i);
  }
  static inline std::vector<int> unmap(const int i, const jackknifeDistribution<T> &from){ 
    return std::vector<int>({i});
  }
};

#endif
