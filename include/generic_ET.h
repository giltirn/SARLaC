#ifndef _GENERIC_ET_H
#define _GENERIC_ET_H

#include<utils.h>

//An expression template engine for any class with an operator()(int) and size() method
template<class T, class Fallback = void>
struct has_ET_tag{ enum{ value = 0 }; };

template<class T>
struct has_ET_tag<T, typename Void<typename T::ET_tag>::type >{ enum{ value = 1 }; };

struct ETleafTag;

template<class T, class Fallback = void>
struct is_ET_leaf{ enum{ value = 0 }; };

template<class T>
struct is_ET_leaf<T, typename Void<typename T::ET_leaf_mark>::type >{ enum{ value = 1 }; };

template<typename T, typename std::enable_if<has_ET_tag<T>::value, int>::type = 0 >
struct get_ET_tag{
  typedef typename std::decay<T>::type::ET_tag type;
};

#define ENABLE_IF_ET_LEAF(T,U) typename std::enable_if<is_ET_leaf<T>::value, U>::type
#define ENABLE_IF_NOT_ET_LEAF(T,U) typename std::enable_if<!is_ET_leaf<T>::value, U>::type
#define ENABLE_IF_TWO_ET_LEAF(T,U,V) typename std::enable_if<is_ET_leaf<T>::value && is_ET_leaf<U>::value, V>::type
#define ENABLE_IF_TWO_ET_LEAF_EQUAL_TAG(T,U,V) typename std::enable_if<is_ET_leaf<T>::value && is_ET_leaf<U>::value && std::is_same<typename T::ET_tag,typename U::ET_tag>::value, V>::type

template<typename T>
struct getElem{
  static inline auto elem(const T &v, const int i)->decltype(v[i]){ return v[i]; }
  static inline auto elem(T &v, const int i)->decltype(v[i]){ return v[i]; }
  static inline size_t common_properties(const T &v){ return v.size(); }
};

//Store references to lvalue operands
template<typename A>
struct ETeval{
  typedef ETleafTag ET_leaf_mark;
  const A &rf;
  typedef ENABLE_IF_NOT_ET_LEAF(A, typename A::ET_tag) ET_tag;
  
  ETeval(const A &r): rf(r){}
  inline decltype(getElem<A>::elem(rf,0)) operator[](const int i) const{ return getElem<A>::elem(rf,i); }
  inline decltype(getElem<A>::common_properties(rf)) common_properties() const{ return getElem<A>::common_properties(rf); }
  static inline auto elem(A &v, const int i)->decltype(getElem<A>::elem(v,i)){ return getElem<A>::elem(v,i); } //used for accessing element in final loop
};

//Store rvalue operands
template<typename A>
struct ETstore{
  typedef ETleafTag ET_leaf_mark;
  A rf;
  typedef ENABLE_IF_NOT_ET_LEAF(A, typename A::ET_tag) ET_tag;
  
  ETstore(A &&r): rf(std::move(r)){  }
  inline decltype(getElem<A>::elem(const_cast<const A &>(rf),0)) operator[](const int i) const{ return getElem<A>::elem(const_cast<const A &>(rf),i); }
  inline decltype(getElem<A>::common_properties(rf)) common_properties() const{ return getElem<A>::common_properties(rf); }
  static inline auto elem(A &v, const int i)->decltype(getElem<A>::elem(v,i)){ return getElem<A>::elem(v,i); } //used for accessing element in final loop
};


template<typename T, typename Fallback = void>
struct rvalueStoreType{
  typedef ETstore<T> type;
};
template<typename T>
struct rvalueStoreType<T, typename Void<typename T::ET_leaf_mark>::type>{
  typedef T type;
};

//Binary operations with both operands vector type
template<template<typename,typename> class Op, typename T, typename U>
struct binaryHelper{
  typedef typename rvalueStoreType<T>::type rT;
  typedef typename rvalueStoreType<U>::type rU;
  
  static inline Op<rT,rU> doit(T &&a, U &&b){ //could be rvalue of operand or of leaf
    return Op<rT,rU>(std::move(a),std::move(b));
  }
  static inline Op<ETeval<T>, ETeval<U> > doit(const T &a, const U &b){
    return Op<ETeval<T>, ETeval<U> >(a,b);
  }
  static inline Op<rT, ETeval<U> > doit(T &&a, const U &b){
    return Op<rT, ETeval<U> >(std::move(a),b);
  }
  static inline Op<ETeval<T>, rU> doit(const T &a, U &&b){
    return Op<ETeval<T>, rU>(a,std::move(b));
  }
};

template<template<typename,typename> class Op, typename Tag>
struct disableGenericETbinOp{
  enum {value = 0};
};

#define ET_BINOP(NAME, OP, OPNAME)		 \
template<typename Leaf1, typename Leaf2> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf1 a; \
  Leaf2 b; \
  typedef ENABLE_IF_TWO_ET_LEAF_EQUAL_TAG(Leaf1,Leaf2, typename Leaf1::ET_tag) ET_tag; \
  \
  NAME(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){ \
    if(a.common_properties() != b.common_properties()) error_exit(std::cout << "Error: binop operands have different common_properties: " << a.common_properties() << ":" << b.common_properties() << std::endl); \
  } \
  \
  inline decltype(a[0] OP b[0]) operator[](const int i) const{ return a[i] OP b[i]; } \
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); } \
};						 \
template<typename T,typename U, \
         typename std::enable_if< \
				  has_ET_tag<typename std::decay<T>::type>::value && \
				  has_ET_tag<typename std::decay<U>::type>::value && \
				  !disableGenericETbinOp<NAME, typename std::decay<T>::type::ET_tag>::value \
				    , int>::type = 0>			\
inline auto OPNAME(T &&a, U &&b)->decltype( binaryHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)) ) { \
  return binaryHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)); \
}

ET_BINOP(ETplus, +, operator+);
ET_BINOP(ETminus, -, operator-);
ET_BINOP(ETtimes, *, operator*);
ET_BINOP(ETdivide, /, operator/);

template<typename T>
struct is_scalar_value{
  enum{ value = std::is_arithmetic<T>::value };
};
template<typename T>
struct is_scalar_value<std::complex<T> >{
  enum{ value = std::is_arithmetic<T>::value };
};

template<typename A>
struct ETscalarEval{
  typedef ETleafTag ET_leaf_mark;
  const A &rf;
  
  ETscalarEval(const A &r): rf(r){}
  inline const A & value() const{ return rf; }
};



template<template<typename,typename> class Op, typename T, typename U>
struct binaryScalarRightHelper{
  typedef typename rvalueStoreType<T>::type rT;
  typedef typename rvalueStoreType<U>::type rU;
  
  static inline Op<rT,ETscalarEval<U> > doit(T &&a, const U &b){
    return Op<rT,ETscalarEval<U> >(std::move(a),b);
  }
  static inline Op<ETeval<T>,ETscalarEval<U> > doit(const T &a, const U &b){
    return Op<ETeval<T>,ETscalarEval<U> >(a,b);
  }
};

#define ET_BINOP_SCALAR_RIGHT(NAME, OP, OPNAME)		 \
template<typename Leaf1, typename Leaf2> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf1 a; \
  Leaf2 b; \
  typedef ENABLE_IF_TWO_ET_LEAF(Leaf1,Leaf2,typename Leaf1::ET_tag) ET_tag; \
  \
  NAME(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){} \
  \
  inline decltype(a[0] OP b.value()) operator[](const int i) const{ return a[i] OP b.value(); } \
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); } \
};   \
template<typename T,typename U, typename std::enable_if<has_ET_tag<typename std::decay<T>::type>::value && is_scalar_value<typename std::decay<U>::type>::value, int>::type = 0> \
inline auto OPNAME(T &&a, U &&b)->decltype(binaryScalarRightHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T >(a),std::forward<U>(b))){ \
  return binaryScalarRightHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)); \
}

ET_BINOP_SCALAR_RIGHT(ETscalarPlusRight, +, operator+);
ET_BINOP_SCALAR_RIGHT(ETscalarMinusRight, -, operator-);
ET_BINOP_SCALAR_RIGHT(ETscalarTimesRight, *, operator*);
ET_BINOP_SCALAR_RIGHT(ETscalarDivideRight, /, operator/);




template<template<typename,typename> class Op, typename T, typename U>
struct binaryScalarLeftHelper{
  typedef typename rvalueStoreType<T>::type rT;
  typedef typename rvalueStoreType<U>::type rU;
  
  static inline Op<ETscalarEval<T>,rU> doit(const T &a, U &&b){
    return Op<ETscalarEval<T>,rU>(a,std::move(b));
  }
  static inline Op<ETscalarEval<T>,ETeval<U> > doit(const T &a, const U &b){
    return Op<ETscalarEval<T>,ETeval<U> >(a,b);
  }
};

#define ET_BINOP_SCALAR_LEFT(NAME, OP, OPNAME)		 \
template<typename Leaf1, typename Leaf2> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf1 a; \
  Leaf2 b; \
  typedef ENABLE_IF_TWO_ET_LEAF(Leaf1,Leaf2,typename Leaf2::ET_tag) ET_tag; \
  \
  NAME(Leaf1 &&aa, Leaf2 &&bb): a(std::move(aa)), b(std::move(bb)){} \
  \
  inline decltype(a.value() OP b[0]) operator[](const int i) const{ return a.value() OP b[i]; } \
  inline decltype(b.common_properties()) common_properties() const{ return b.common_properties(); } \
};   \
template<typename T,typename U, typename std::enable_if<is_scalar_value<typename std::decay<T>::type>::value && has_ET_tag<typename std::decay<U>::type>::value , int>::type = 0> \
inline auto OPNAME(T &&a, U &&b)->decltype(binaryScalarLeftHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T >(a),std::forward<U>(b))){ \
  return binaryScalarLeftHelper<NAME,typename std::decay<T>::type,typename std::decay<U>::type>::doit(std::forward<T>(a),std::forward<U>(b)); \
}

ET_BINOP_SCALAR_LEFT(ETscalarPlusLeft, +, operator+);
ET_BINOP_SCALAR_LEFT(ETscalarMinusLeft, -, operator-);
ET_BINOP_SCALAR_LEFT(ETscalarTimesLeft, *, operator*);




template<template<typename> class Op, typename T>
struct unaryHelper{
  typedef typename rvalueStoreType<T>::type rT;
  
  static inline Op<ETeval<T> > doit(const T &a){
    return Op<ETeval<T> >(a);
  }
  static inline Op<rT> doit(T &&a){
    return Op<rT>(std::move(a));
  }
};
#define ET_UNOP(NAME, OP, OPNAME)		 \
template<typename Leaf> \
struct NAME{ \
  typedef ETleafTag ET_leaf_mark; \
  Leaf a; \
  typedef ENABLE_IF_ET_LEAF(Leaf, typename Leaf::ET_tag) ET_tag;	\
  \
  NAME(Leaf &&aa): a(std::move(aa)){} \
  \
  inline decltype(OP(a[0])) operator[](const int i) const{ return OP(a[i]); } \
  inline decltype(a.common_properties()) common_properties() const{ return a.common_properties(); } \
};   \
template<typename T, typename std::enable_if<has_ET_tag<typename std::decay<T>::type>::value , int>::type = 0> \
inline auto OPNAME(T &&a)->decltype(unaryHelper<NAME,typename std::decay<T>::type>::doit(std::forward<T >(a))){ \
  return unaryHelper<NAME,typename std::decay<T>::type>::doit(std::forward<T>(a)); \
}

template<typename T>
inline T negate(const T &a){ return -a; }

ET_UNOP(ETnegate, negate, operator-);
ET_UNOP(ETexp, exp, exp);
ET_UNOP(ETsqrt, sqrt, sqrt);

//PARALLELIZE_DISTRIBUTION_ET disabled by default because it makes performance worse in most cases (better with very larger, O(10000) configs)
#ifdef PARALLELIZE_DISTRIBUTION_ET
template<typename T, typename U>
inline void loop_eval(T &obj, U &&expr){
#pragma omp parallel for
    for(int i=0;i<obj.size();i++) ETeval<T>::elem(obj,i) = expr[i];
}
#else
template<typename T, typename U>
inline void loop_eval(T &obj, U &&expr){
    for(int i=0;i<obj.size();i++) ETeval<T>::elem(obj,i) = expr[i];
}
#endif

//Put this inside your class to enable the ET
//Tag is used to discriminate between classes of object; a binary op requires both ops have the same tag
#define ENABLE_GENERIC_ET(CLASS_NAME, CLASS_FULL, TAG)			\
  typedef TAG ET_tag;							\
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,CLASS_FULL>::value, int>::type = 0> \
  CLASS_NAME(U&& expr): CLASS_NAME(expr.common_properties()){			\
    loop_eval<CLASS_FULL, U>(*this, std::forward<U>(expr));		\
  }\
  \
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,CLASS_FULL>::value, int>::type = 0> \
  CLASS_FULL & operator=(U && expr){ \
    this->resize(expr.common_properties());    \
    loop_eval<CLASS_FULL, U>(*this, std::forward<U>(expr));	\
    return *this; \
  }



#endif
