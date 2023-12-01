#ifndef _SARLAC_RANDOMNUMBER_H
#define _SARLAC_RANDOMNUMBER_H

//Generate random real and complex numbers from different distributions

#include <complex>

#include<config.h>
#include<utils/template_wizardry.h>
#include<random/rng.h>

SARLAC_START_NAMESPACE

//Uniform random numbers
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type uniformRandom(T &v, const T start, const T end, RNGstore &rng = RNG){
  std::uniform_real_distribution<> dis(start, end);
  v = dis(rng());
}
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type uniformRandom(std::complex<T> &v, const T start, const T end, RNGstore &rng = RNG){
  std::uniform_real_distribution<> dis(start, end);
  reinterpret_cast<T(&)[2]>(v)[0] = dis(rng());
  reinterpret_cast<T(&)[2]>(v)[1] = dis(rng());
}
template<typename T>
inline typename std::enable_if<hasSampleMethod<T>::value && hasDataType<T>::value, void>::type uniformRandom(T &v, const typename T::DataType start, const typename T::DataType end, RNGstore &rng = RNG){ //for distribution types
  for(int i=0;i<v.size();i++) uniformRandom(v.sample(i),start,end,rng);
}

//Gaussian random numbers
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type gaussianRandom(T &v, const T mean, const T stddev, RNGstore &rng = RNG){
  std::normal_distribution<> dis(mean, stddev);
  v = dis(rng());
}
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type gaussianRandom(std::complex<T> &v, const T mean, const T stddev, RNGstore &rng = RNG){ //real and complex from same distribution. 
  std::normal_distribution<> dis(mean, stddev);
  reinterpret_cast<T(&)[2]>(v)[0] = dis(rng());
  reinterpret_cast<T(&)[2]>(v)[1] = dis(rng());
}
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type gaussianRandom(std::complex<T> &v, const std::complex<T> mean, const std::complex<T> stddev, RNGstore &rng = RNG){ //real and complex from same distribution. 
  std::normal_distribution<> redis(mean.real(), stddev.real());
  std::normal_distribution<> imdis(mean.imag(), stddev.imag());
  reinterpret_cast<T(&)[2]>(v)[0] = redis(rng());
  reinterpret_cast<T(&)[2]>(v)[1] = redis(rng());
}


template<typename T>
inline typename std::enable_if<hasSampleMethod<T>::value && hasDataType<T>::value, void>::type gaussianRandom(T &v, const typename T::DataType mean, const typename T::DataType stddev, RNGstore &rng = RNG){ //for distribution types
  for(int i=0;i<v.size();i++) gaussianRandom(v.sample(i),mean,stddev,rng);
}
template<typename T>
inline void gaussianRandom(std::vector<T> &v, const T mean, const T stddev, RNGstore &rng = RNG){
  for(int i=0;i<v.size();i++) gaussianRandom(v[i],mean,stddev,rng);
}
template<typename T>
inline void gaussianRandom(std::vector<std::complex<T> > &v, const std::complex<T> mean, const std::complex<T> stddev, RNGstore &rng = RNG){
  for(int i=0;i<v.size();i++) gaussianRandom(v[i],mean,stddev,rng);
}
template<typename T>
inline void gaussianRandom(std::vector<std::complex<T> > &v, const T mean, const T stddev, RNGstore &rng = RNG){
  for(int i=0;i<v.size();i++) gaussianRandom(v[i],mean,stddev,rng);
}


//Alternative calling convention that returns the random number
template<typename T>
struct _rangetype{ typedef T type; };
template<typename T>
struct _rangetype<std::complex<T> >{ typedef T type; };


template<typename T>
inline T uniformRandom(const typename _rangetype<T>::type start, const typename _rangetype<T>::type end, RNGstore &rng = RNG){
  T out;
  uniformRandom(out,start,end,rng);
  return out;
}
template<typename T>
inline T gaussianRandom(const typename _rangetype<T>::type mean, const typename _rangetype<T>::type stddev, RNGstore &rng = RNG){
  T out;
  gaussianRandom(out,mean,stddev,rng);
  return out;
}

SARLAC_END_NAMESPACE

#endif
