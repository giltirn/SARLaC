#ifndef _CPSFIT_UTILS_MATH_H_
#define _CPSFIT_UTILS_MATH_H_

//Some math-related functionality
#include<complex>
#include<omp.h>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry/distribution.h>
#include<utils/template_wizardry/other.h>
#include<utils/template_wizardry/matrix_vector.h>

CPSFIT_START_NAMESPACE

//Extract real and imaginary parts of complex by reference
template<typename T>
inline T & real_ref(std::complex<T> &z){
  return reinterpret_cast<T(&)[2]>(z)[0];
}
template<typename T>
inline T & imag_ref(std::complex<T> &z){
  return reinterpret_cast<T(&)[2]>(z)[1];
}


//Set a type to zero. Works for types with assignment to double method, zero method or basic POD floating point types.
template<typename T, typename std::enable_if<hasZeroMethod<T>::value, int>::type = 0>
inline void zeroit(T &v){
  v.zero();
}
template<typename T, typename std::enable_if<!hasZeroMethod<T>::value && (hasEqualsMethod<T,double>::value || hasEqualsMethod<T,float>::value), int>::type = 0>
inline void zeroit(T &v){
  v = 0.;
}
inline void zeroit(double &v){
  v = 0.;
}
inline void zeroit(float &v){
  v = 0.;
}

///Perform the threaded sum of an operator struct with operator()(const int)  and size() methods
template<typename Operation, typename std::enable_if<hasParenthesesConstAccessor<Operation>::value,int>::type = 0>
auto threadedSum(const Operation &op)->typename std::decay<decltype(op(0))>::type{
  typedef typename std::decay<decltype(op(0))>::type T;
  const int N = op.size();
  const int nthread = omp_get_max_threads();
 
  T init_zero(op(0)); zeroit(init_zero);    //_threadedSumHelper<Operation,T>::getZero(op);
  std::vector<T> sum(nthread, init_zero);

#pragma omp parallel for
  for(int i=0;i<N;i++)
    sum[omp_get_thread_num()] = sum[omp_get_thread_num()] + op(i);

  for(int t=1;t<nthread;t++)
    sum[0] = sum[0] + sum[t];

  return sum[0];
}
					
//Perform the threaded sum of the elements of a vector-like array (operator[](const int) and size() methods)	    
template<typename VectorOfData, typename std::enable_if<hasSquareBracketsConstAccessor<VectorOfData>::value,int>::type = 0>
typename std::decay<decltype(VectorOfData()[0])>::type threadedSum(const VectorOfData &v){
  typedef typename std::decay<decltype(VectorOfData()[0])>::type T;
  struct Op{
    const VectorOfData &vv;
    inline const T & operator()(const int idx) const{ return vv[idx]; }
    inline size_t size() const{ return vv.size(); }
    Op(const VectorOfData &_vv): vv(_vv){}
  };
  Op op(v);
  return threadedSum<Op>(op);
}


CPSFIT_END_NAMESPACE
#endif
