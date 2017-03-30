#ifndef _NUMERIC_TENSORS_ET_H_
#define _NUMERIC_TENSORS_ET_H_

#include<template_wizardry.h>

template<typename T>
struct VectorType{};
template<typename T>
struct MatrixType{};

template<class T>
struct is_VectorType{ enum{value = 0}; };

template<class T>
struct is_VectorType< VectorType<T> >{ enum{value = 1}; };

template<class T>
struct is_MatrixType{ enum{value = 0}; };

template<class T>
struct is_MatrixType< MatrixType<T> >{ enum{value = 1}; };


template<class T, class Fallback = void>
struct base_type_is_vector{ enum{ value = 0 }; };

template<class T>
struct base_type_is_vector<T, typename Void<typename T::Tensor_ET_base_type>::type >{ enum{ value = is_VectorType<typename T::Tensor_ET_base_type>::value }; };

template<class T, class Fallback = void>
struct base_type_is_matrix{ enum{ value = 0 }; };

template<class T>
struct base_type_is_matrix<T, typename Void<typename T::Tensor_ET_base_type>::type >{ enum{ value = is_MatrixType<typename T::Tensor_ET_base_type>::value }; };

#define IS_EXPRESSION_WITH_VECTOR_BASE_TYPE(EXPR, ME) typename std::enable_if<  \
  !std::is_same<EXPR,ME>::value && \
  base_type_is_vector<EXPR>::value	   \
  , int>::type = 0

#define IS_EXPRESSION_WITH_MATRIX_BASE_TYPE(EXPR, ME) typename std::enable_if<  \
  !std::is_same<EXPR,ME>::value && \
  base_type_is_matrix<EXPR>::value	   \
  , int>::type = 0

//Vector binary operators
#define TYPES_VECTOR_VECTOR(T,U) typename std::enable_if<base_type_is_vector<T>::value && base_type_is_vector<U>::value, int>::type = 0
#define TYPES_VECTOR_SCALAR(T,U) typename std::enable_if<base_type_is_vector<T>::value && is_scalar<U>::value, int>::type = 0

#define DEF_VECTOR_BINOP_VV(CLASS, OP, OP_ACTION_AB) \
template<typename T, typename U, TYPES_VECTOR_VECTOR(T,U) > \
class CLASS{ \
public: \
  typedef typename T::Tensor_ET_base_type Tensor_ET_base_type; \
  T const& a; \
  U const& b; \
  CLASS(const T &aa, const U &bb): a(aa), b(bb){ assert(aa.size() == bb.size()); } \
  \
  inline auto operator[](const int i) const ->decltype(OP_ACTION_AB){ return OP_ACTION_AB; } \
  inline int size() const{ return a.size(); } \
}; \
\
template<typename T, typename U, TYPES_VECTOR_VECTOR(T,U) >  \
CLASS<T,U> OP(const T &a, const U &b){		\
  return CLASS<T,U>(a,b); \
}  

DEF_VECTOR_BINOP_VV(ETvectorPlus, operator+, a[i] + b[i]);
DEF_VECTOR_BINOP_VV(ETvectorMinus, operator-, a[i] - b[i]);
DEF_VECTOR_BINOP_VV(ETvectorTimes, operator*, a[i] * b[i]);
DEF_VECTOR_BINOP_VV(ETvectorDivide, operator/, a[i] / b[i]);

#define DEF_VECTOR_BINOP_VS(CLASS, OP, OP_ACTION_AB) \
template<typename T, typename U, TYPES_VECTOR_SCALAR(T,U) > \
class CLASS{ \
public: \
  typedef typename T::Tensor_ET_base_type Tensor_ET_base_type; \
  T const& a; \
  U const& b; \
  CLASS(const T &aa, const U &bb): a(aa), b(bb){ } \
  \
  inline auto operator[](const int i) const ->decltype(OP_ACTION_AB){ return OP_ACTION_AB; } \
  inline int size() const{ return a.size(); } \
};

#define DEF_VECTOR_BINOP_VS_VOPS(CLASS,OP)		     \
template<typename T, typename U, TYPES_VECTOR_SCALAR(T,U) >  \
CLASS<T,U> OP(const T &a, const U &b){		\
  return CLASS<T,U>(a,b); \
}

#define DEF_VECTOR_BINOP_VS_SOPV(CLASS,OP)		     \
template<typename T, typename U, TYPES_VECTOR_SCALAR(U,T) >  \
CLASS<U,T> OP(const T &a, const U &b){		\
  return CLASS<U,T>(b,a); \
}

DEF_VECTOR_BINOP_VS(ETvectorScalarMult, operator*, a[i]*b);
DEF_VECTOR_BINOP_VS_VOPS(ETvectorScalarMult, operator*);
DEF_VECTOR_BINOP_VS_SOPV(ETvectorScalarMult, operator*);

DEF_VECTOR_BINOP_VS(ETvectorScalarDivide, operator/, a[i]/b);
DEF_VECTOR_BINOP_VS_VOPS(ETvectorScalarDivide, operator/);


//Matrix binary operators
#define TYPES_MATRIX_MATRIX(T,U) typename std::enable_if<base_type_is_matrix<T>::value && base_type_is_matrix<U>::value, int>::type = 0
#define TYPES_MATRIX_SCALAR(T,U) typename std::enable_if<base_type_is_matrix<T>::value && is_scalar<U>::value, int>::type = 0



//Use function model for matrix because the operations aren't as uniform as vector
#define DEF_MATRIX_BINOP_MM(CLASS, OP, OP_ACTION_AB_FUNC) \
template<typename T, typename U, TYPES_MATRIX_MATRIX(T,U) > \
class CLASS{ \
public: \
  typedef typename T::Tensor_ET_base_type Tensor_ET_base_type; \
  T const& a; \
  U const& b; \
  CLASS(const T &aa, const U &bb): a(aa), b(bb){ assert(aa.size() == bb.size()); } \
  \
  inline auto operator()(const int i, const int j) const ->decltype(OP_ACTION_AB_FUNC(i,j,a,b)){ return OP_ACTION_AB_FUNC(i,j,a,b); } \
  inline int size() const{ return a.size(); } \
}; \
\
template<typename T, typename U, TYPES_MATRIX_MATRIX(T,U) >  \
CLASS<T,U> OP(const T &a, const U &b){		\
  return CLASS<T,U>(a,b); \
}  

template<typename T,typename U>
inline auto _mat_mat_sum(const int i, const int j, const T &a, const U &b)->decltype(a(0,0)+a(0,0)){ return a(i,j)+b(i,j); }
DEF_MATRIX_BINOP_MM(ETmatrixPlus, operator+, _mat_mat_sum);

#endif
