#ifndef _CPSFIT_RANDOM_H
#define _CPSFIT_RANDOM_H

#include <random>
#include <complex>

class RNGstore{
public:
  typedef std::mt19937 RNGtype;
  typedef typename RNGtype::result_type seedType;
private:
  RNGtype* rng;
public:
  RNGstore(): rng(NULL){}

  void initialize(const seedType seed){ if(rng==NULL) rng = new RNGtype(seed); }
  void initialize(){ if(rng == NULL) rng = new RNGtype(); }
  
  RNGstore(const seedType seed){ initialize(seed); }

  bool isInitialized() const{ return rng != NULL; }
  
  RNGtype &operator()(){ return *rng; }
};

RNGstore RNG; //static instance

template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type uniformRandom(T &v, const T start, const T end){
  std::uniform_real_distribution<> dis(start, end);
  v = dis(RNG());
}
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type uniformRandom(std::complex<T> &v, const T start, const T end){
  std::uniform_real_distribution<> dis(start, end);
  reinterpret_cast<T(&)[2]>(v)[0] = dis(RNG());
  reinterpret_cast<T(&)[2]>(v)[1] = dis(RNG());
}
template<typename T>
inline typename std::enable_if<hasSampleMethod<T>::value && hasDataType<T>::value, void>::type uniformRandom(T &v, const typename T::DataType start, const typename T::DataType end){ //for distribution types
  for(int i=0;i<v.size();i++) uniformRandom(v.sample(i),start,end);
}

template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type gaussianRandom(T &v, const T mean, const T stddev){
  std::normal_distribution<> dis(mean, stddev);
  v = dis(RNG());
}
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type gaussianRandom(std::complex<T> &v, const T mean, const T stddev){ //real and complex from same distribution. 
  std::normal_distribution<> dis(mean, stddev);
  reinterpret_cast<T(&)[2]>(v)[0] = dis(RNG());
  reinterpret_cast<T(&)[2]>(v)[1] = dis(RNG());
}
template<typename T>
inline typename std::enable_if<hasSampleMethod<T>::value && hasDataType<T>::value, void>::type gaussianRandom(T &v, const typename T::DataType mean, const typename T::DataType stddev){ //for distribution types
  for(int i=0;i<v.size();i++) gaussianRandom(v.sample(i),mean,stddev);
}


template<typename T>
struct _rangetype{ typedef T type; };
template<typename T>
struct _rangetype<std::complex<T> >{ typedef T type; };


template<typename T>
inline T uniformRandom(const typename _rangetype<T>::type start, const typename _rangetype<T>::type end){
  T out;
  uniformRandom(out,start,end);
  return out;
}
template<typename T>
inline T gaussianRandom(const typename _rangetype<T>::type mean, const typename _rangetype<T>::type stddev){
  T out;
  gaussianRandom(out,mean,stddev);
  return out;
}

#endif
