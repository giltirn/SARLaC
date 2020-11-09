#ifndef _CPSFIT_NUMERIC_VECTOR_MATH_H_
#define _CPSFIT_NUMERIC_VECTOR_MATH_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector/class.h>
#include<utils/template_wizardry/complexify.h>

CPSFIT_START_NAMESPACE

//Vector dot product
template<typename T>
T dot(const NumericVector<T> &a, const NumericVector<T> &b){
  assert(a.size() == b.size());
  T out(a(0)); zeroit(out);
  for(int i=0;i<a.size();i++) out = out + a(i) * b(i);
  return out;
}

//Complex dot product
template<typename T>
T complex_dot(const NumericVector<T> &a, const NumericVector<T> &b){
  assert(a.size() == b.size());
  T out(a(0)); zeroit(out);
  for(int i=0;i<a.size();i++) out = out + conj(a(i)) * b(i);
  return out;
}

//Vector modulus
template<typename T>
T mod2(const NumericVector<T> &m){
  return dot(m,m);
}

//Complex modulus
template<typename T>
Realify<T> complex_mod2(const NumericVector<T> &m){
  return real(complex_dot(m,m));
}

//Generate an orthonormal basis from a set of vectors
template<typename T>
std::vector<NumericVector<T> > GrammSchmidtOrthonormalize(const std::vector<NumericVector<T> > &v){
  int n = v.size();
  int L = v[0].size();
  for(int i=1;i<n;i++) assert(v[i].size() == L);

  std::vector<NumericVector<T> > u(n);
  u[0] = v[0]/sqrt(mod2(v[0]));

  for(int i=1;i<n;i++){
    u[i] = v[i];
    for(int j=0; j<i; j++)
      u[i] = u[i] - dot(v[i], u[j])*u[j];
    u[i] = u[i] / sqrt(mod2(u[i]));
  }
  return u;
}


CPSFIT_END_NAMESPACE

#endif
