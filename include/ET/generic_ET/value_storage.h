#ifndef _GENERIC_ET_VALUE_STORAGE_H_
#define _GENERIC_ET_VALUE_STORAGE_H_

#include<ET/generic_ET/tagging.h>
#include<ET/generic_ET/getelem.h>

CPSFIT_START_NAMESPACE

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



template<typename T>
struct is_scalar_value{
  enum{ value = std::is_arithmetic<T>::value };
};
template<typename T>
struct is_scalar_value<std::complex<T> >{
  enum{ value = std::is_arithmetic<T>::value };
};

template<typename A>
struct ETscalarEval{ //store scalar references
  typedef ETleafTag ET_leaf_mark;
  const A &rf;
  
  ETscalarEval(const A &r): rf(r){}
  inline const A & value() const{ return rf; }
};

CPSFIT_END_NAMESPACE

#endif
