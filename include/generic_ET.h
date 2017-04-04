#ifndef _GENERIC_ET_H
#define _GENERIC_ET_H

//An expression template engine for any class with an operator()(int) and size() method


template<class T, class BaseType, class Fallback = void>
struct generic_ET_base_type_equals{ enum{ value = 0 }; };

template<class T, class BaseType>
struct generic_ET_base_type_equals<T, BaseType, typename Void<typename T::Generic_ET_base_type>::type >{ enum{ value = std::is_same<typename T::Generic_ET_base_type, BaseType>::value }; };

#define BASE_TYPES_EQUAL_TO(T,U,BaseType) typename std::enable_if<generic_ET_base_type_equals<T,BaseType>::value && generic_ET_base_type_equals<U,BaseType>::value, int>::type = 0

#define BASE_TYPE_EQUAL_TO_VECTORSCALAR(T,U,BaseType) typename std::enable_if<generic_ET_base_type_equals<T,BaseType>::value && is_scalar<U>::value, int>::type = 0

#define BASE_TYPE_EQUAL_TO_UNOP(T,BaseType) typename std::enable_if<generic_ET_base_type_equals<T,BaseType>::value, int>::type = 0

#define IS_EXPRESSION_WITH_GENERIC_BASE_TYPE(EXPR, ME, BASETYPE) typename std::enable_if< \
  !std::is_same<EXPR,ME>::value && \
  generic_ET_base_type_equals<EXPR, BASETYPE>::value	\
  , int>::type = 0


#define DEF_GENERIC_BINOP(CLASS, OP, OP_ACTION_AB, BASETYPE)	    \
  template<typename T, typename U, BASE_TYPES_EQUAL_TO(T,U,BASETYPE) >	    \
class CLASS{ \
public: \
  typedef typename T::Generic_ET_base_type Generic_ET_base_type; \
  T const& a; \
  U const& b; \
  CLASS(const T &aa, const U &bb): a(aa), b(bb){ assert(aa.size() == bb.size()); } \
  \
  inline auto operator()(const int i) const ->decltype(OP_ACTION_AB){ return OP_ACTION_AB; } \
  inline int size() const{ return a.size(); } \
}; \
\
template<typename T, typename U, BASE_TYPES_EQUAL_TO(T,U,BASETYPE) >  \
CLASS<T,U> OP(const T &a, const U &b){		\
  return CLASS<T,U>(a,b); \
}  

//The below are broken into a base and 2 ops because it doesn't always make sense to perform a vector-scalar binop in both orders  eg  vector/scalar but not scalar/vector

//vector-scalar type
#define DEF_GENERIC_BINOP_VS_BASE(CLASS, OP, OP_ACTION_AB, BASETYPE) \
template<typename T, typename U, BASE_TYPE_EQUAL_TO_VECTORSCALAR(T,U,BASETYPE) > \
class CLASS{ \
public: \
  typedef typename T::Generic_ET_base_type Generic_ET_base_type; \
  T const& a; \
  U const& b; \
  CLASS(const T &aa, const U &bb): a(aa), b(bb){ } \
  \
  inline auto operator()(const int i) const ->decltype(OP_ACTION_AB){ return OP_ACTION_AB; } \
  inline int size() const{ return a.size(); } \
};

//T is vector-type, U is scalar type
#define DEF_GENERIC_BINOP_VS(CLASS,OP, BASETYPE)		     \
template<typename T, typename U, BASE_TYPE_EQUAL_TO_VECTORSCALAR(T,U,BASETYPE) >  \
CLASS<T,U> OP(const T &a, const U &b){		\
  return CLASS<T,U>(a,b); \
}

//T is scalar-type, U is vector type
#define DEF_GENERIC_BINOP_SV(CLASS,OP, BASETYPE)		     \
template<typename T, typename U, BASE_TYPE_EQUAL_TO_VECTORSCALAR(U,T,BASETYPE) >  \
CLASS<U,T> OP(const T &a, const U &b){		\
  return CLASS<U,T>(b,a); \
}


#define DEF_GENERIC_UNOP(CLASS, OP, OP_ACTION, BASETYPE) \
template<typename T, BASE_TYPE_EQUAL_TO_UNOP(T,BASETYPE) > \
class CLASS{ \
public: \
  typedef typename T::Generic_ET_base_type Generic_ET_base_type; \
  T const& a; \
  CLASS(const T &aa): a(aa){ } \
  \
  inline auto operator()(const int i) const ->decltype(OP_ACTION){ return OP_ACTION; } \
  inline int size() const{ return a.size(); } \
}; \
template<typename T, BASE_TYPE_EQUAL_TO_UNOP(T,BASETYPE) >  \
CLASS<T> OP(const T &a){		\
  return CLASS<T>(a); \
}

#define PREPROC_CONCAT_(A,B) A ## B
#define PREPROC_CONCAT(A,B) PREPROC_CONCAT_(A,B)

//Put this inside your class. TAG is a struct type used to differentiate the ETE for this type from others. Most often convenient to set equal to the class name
#define ENABLE_GENERIC_ET(CLASS, TAG) \
  typedef Params Generic_ET_base_type; \
  \
  template<typename Expr, IS_EXPRESSION_WITH_GENERIC_BASE_TYPE(Expr, CLASS, TAG) > \
  CLASS(const Expr &e){ \
    for(int i=0;i<this->size();i++) this->operator()(i) = e(i); \
  }

//Put this outside your struct
#define DEF_GENERIC_ET(TAG) \
  DEF_GENERIC_BINOP( PREPROC_CONCAT(ETplus_,TAG), operator+, a(i) + b(i), TAG ); \
  DEF_GENERIC_BINOP( PREPROC_CONCAT(ETminus_,TAG), operator-, a(i) - b(i), TAG ); \
  DEF_GENERIC_BINOP( PREPROC_CONCAT(ETtimes_,TAG), operator*, a(i) * b(i), TAG ); \
  DEF_GENERIC_BINOP( PREPROC_CONCAT(ETdivide_,TAG), operator/, a(i) / b(i), TAG );	\
  \
  DEF_GENERIC_BINOP_VS_BASE( PREPROC_CONCAT(ETscalarMult_,TAG), operator*, a(i)*b, TAG); \
  DEF_GENERIC_BINOP_VS( PREPROC_CONCAT(ETscalarMult_,TAG), operator*, TAG); \
  DEF_GENERIC_BINOP_SV( PREPROC_CONCAT(ETscalarMult_,TAG), operator*, TAG); \
  \
  DEF_GENERIC_BINOP_VS_BASE( PREPROC_CONCAT(ETscalarDiv_,TAG), operator/, a(i)/b, TAG); \
  DEF_GENERIC_BINOP_VS( PREPROC_CONCAT(ETscalarDiv_,TAG), operator/, TAG); \
  \
  DEF_GENERIC_UNOP( PREPROC_CONCAT(ETsqrt_,TAG), sqrt, sqrt(a(i)), TAG)




#endif
