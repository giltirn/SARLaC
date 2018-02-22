#ifndef _CPSFIT_RANDOMNUMBER_H
#define _CPSFIT_RANDOMNUMBER_H

//Generate random real and complex numbers from different distributions

#include <complex>

#include<config.h>
#include<utils/template_wizardry.h>
#include<random/rng.h>

CPSFIT_START_NAMESPACE

//Uniform random numbers
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

//Gaussian random numbers
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
inline typename std::enable_if<std::is_floating_point<T>::value, void>::type gaussianRandom(std::complex<T> &v, const std::complex<T> mean, const std::complex<T> stddev){ //real and complex from same distribution. 
  std::normal_distribution<> redis(mean.real(), stddev.real());
  std::normal_distribution<> imdis(mean.imag(), stddev.imag());
  reinterpret_cast<T(&)[2]>(v)[0] = redis(RNG());
  reinterpret_cast<T(&)[2]>(v)[1] = redis(RNG());
}


template<typename T>
inline typename std::enable_if<hasSampleMethod<T>::value && hasDataType<T>::value, void>::type gaussianRandom(T &v, const typename T::DataType mean, const typename T::DataType stddev){ //for distribution types
  for(int i=0;i<v.size();i++) gaussianRandom(v.sample(i),mean,stddev);
}
template<typename T>
inline void gaussianRandom(std::vector<T> &v, const T mean, const T stddev){
  for(int i=0;i<v.size();i++) gaussianRandom(v[i],mean,stddev);
}
template<typename T>
inline void gaussianRandom(std::vector<std::complex<T> > &v, const std::complex<T> mean, const std::complex<T> stddev){
  for(int i=0;i<v.size();i++) gaussianRandom(v[i],mean,stddev);
}
template<typename T>
inline void gaussianRandom(std::vector<std::complex<T> > &v, const T mean, const T stddev){
  for(int i=0;i<v.size();i++) gaussianRandom(v[i],mean,stddev);
}


//Alternative calling convention that returns the random number
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

CPSFIT_END_NAMESPACE

#endif
