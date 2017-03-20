#ifndef _DISTRIBUTION_ET_H
#define _DISTRIBUTION_ET_H

#include<cstdio>
#include<vector>
#include<type_traits>
#include<cassert>
#include<template_wizardry.h>


template<class T, class Fallback = void>
struct has_ET_base_type{ enum{ value = 0 }; };

template<class T>
struct has_ET_base_type<T, typename Void<typename T::ET_base_type>::type >{ enum{ value = 1 }; };

#define HAS_ET_BASE_TYPE(T) typename std::enable_if<has_ET_base_type<T>::value, int>::type = 0
#define SAME_BASE_TYPE(T,U) typename std::enable_if<std::is_same<typename T::ET_base_type, typename U::ET_base_type>::value, int>::type = 0

#define IS_EXPRESSION_WITH_SAME_BASE_TYPE(EXPR, ME) typename std::enable_if<  \
!std::is_same<EXPR,ME>::value &&   \
has_ET_base_type<EXPR>::value &&   \
has_ET_base_type<ME>::value &&     \
std::is_same<typename EXPR::ET_base_type, typename ME::ET_base_type>::value \
  , int>::type = 0


//Binary operators
template<typename T, typename U, SAME_BASE_TYPE(T,U) >
class ETplus{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  U const& b;
  ETplus(const T &aa, const U &bb): a(aa), b(bb){ assert(aa.size() == bb.size()); }
  
  inline decltype(a.sample(0) + b.sample(0)) sample(const int i) const{ return a.sample(i) + b.sample(i); }
  inline int size() const{ return a.size(); }
};

template<typename T, typename U, SAME_BASE_TYPE(T,U) > 
ETplus<T,U> operator+(const T &a, const U &b){
  return ETplus<T,U>(a,b);
}


template<typename T, typename U, SAME_BASE_TYPE(T,U) >
class ETminus{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  U const& b;
  ETminus(const T &aa, const U &bb): a(aa), b(bb){ assert(aa.size() == bb.size()); }
  
  inline decltype(a.sample(0) - b.sample(0)) sample(const int i) const{ return a.sample(i) - b.sample(i); }
  inline int size() const{ return a.size(); }
};

template<typename T, typename U, SAME_BASE_TYPE(T,U) > 
ETminus<T,U> operator-(const T &a, const U &b){
  return ETminus<T,U>(a,b);
}


template<typename T, typename U, SAME_BASE_TYPE(T,U) >
class ETmult{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  U const& b;
  ETmult(const T &aa, const U &bb): a(aa), b(bb){ assert(aa.size() == bb.size()); }
  
  inline decltype(a.sample(0) * b.sample(0)) sample(const int i) const{ return a.sample(i) * b.sample(i); }
  inline int size() const{ return a.size(); }
};

template<typename T, typename U, SAME_BASE_TYPE(T,U) > 
ETmult<T,U> operator*(const T &a, const U &b){
  return ETmult<T,U>(a,b);
}

template<typename T, typename U, SAME_BASE_TYPE(T,U) >
class ETdiv{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  U const& b;
  ETdiv(const T &aa, const U &bb): a(aa), b(bb){ assert(aa.size() == bb.size()); }
  
  inline decltype(a.sample(0) / b.sample(0)) sample(const int i) const{ return a.sample(i) / b.sample(i); }
  inline int size() const{ return a.size(); }
};

template<typename T, typename U, SAME_BASE_TYPE(T,U) > 
ETdiv<T,U> operator/(const T &a, const U &b){
  return ETdiv<T,U>(a,b);
}

template<typename T, typename Scalar>
class ETscalarMult{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  Scalar const& b;
  ETscalarMult(const T &aa, const Scalar &bb): a(aa), b(bb){ }
  
  inline decltype(a.sample(0) * b) sample(const int i) const{ return a.sample(i)*b; }
  inline int size() const{ return a.size(); }
};

template<typename T, HAS_ET_BASE_TYPE(T)>
ETscalarMult<T,double> operator*(const T &a, const double &b){
  return ETscalarMult<T,double>(a,b);
}
template<typename T, HAS_ET_BASE_TYPE(T)>
ETscalarMult<T,double> operator*(const double &a, const T &b){
  return ETscalarMult<T,double>(b,a);
}


template<typename T, typename Scalar>
class ETscalarDiv{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  Scalar const& b;
  ETscalarDiv(const T &aa, const Scalar &bb): a(aa), b(bb){ }
  
  inline decltype(a.sample(0) / b) sample(const int i) const{ return a.sample(i)/b; }
  inline int size() const{ return a.size(); }
};

template<typename T, HAS_ET_BASE_TYPE(T)>
ETscalarDiv<T,double> operator/(const T &a, const double &b){
  return ETscalarDiv<T,double>(a,b);
}


template<typename Scalar,typename T>
class ETscalarNumDiv{
public:
  typedef typename T::ET_base_type ET_base_type;
  Scalar const& a;
  T const& b;
  ETscalarNumDiv(const Scalar &aa, const T &bb): a(aa), b(bb){ }
  
  inline decltype(a / b.sample(0)) sample(const int i) const{ return a/b.sample(i); }
  inline int size() const{ return b.size(); }
};

template<typename T, HAS_ET_BASE_TYPE(T)>
ETscalarNumDiv<T,double> operator/(const double &a, const T &b){
  return ETscalarNumDiv<double,T>(a,b);
}




//Unary operators

template<typename T>
class ETsqrt{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  ETsqrt(const T &aa): a(aa){}
  
  inline decltype(sqrt(a.sample(0)) ) sample(const int i) const{ return sqrt(a.sample(i)); }
  inline int size() const{ return a.size(); }
};

template<typename T,  HAS_ET_BASE_TYPE(T)> 
ETsqrt<T> sqrt(const T &a){
  return ETsqrt<T>(a);
}


template<typename T>
class ETneg{
public:
  typedef typename T::ET_base_type ET_base_type;
  T const& a;
  ETneg(const T &aa): a(aa){}
  
  inline decltype(-a.sample(0)) sample(const int i) const{ return -a.sample(i); }
  inline int size() const{ return a.size(); }
};

template<typename T,  HAS_ET_BASE_TYPE(T)> 
ETneg<T> operator-(const T &a){
  return ETneg<T>(a);
}




#endif
