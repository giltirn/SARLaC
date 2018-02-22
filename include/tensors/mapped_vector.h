#ifndef _CPSFIT_MAPPED_VECTOR_H_
#define _CPSFIT_MAPPED_VECTOR_H_

//A vector type with an additional global mapping between some generic tag and it's elements

#include<map>

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>

CPSFIT_START_NAMESPACE

template<typename T, typename MapType>
class mappedVector: public NumericVector<T>{
  //MapType must conform to the following form
  // struct MapTypeExample{
  //   typedef std::string tagType;    //choose a tag type - here a string but can be anything
  //   inline int map(const tagType &tag) const; //mapping between tag and index
  //   inline tagType unmap(const int idx) const; //unmapping between index and tag
  //   inline int size() const; //number of elements
  // };
  
  //An example is stringTagMap

  MapType const* mapping;

public:
  typedef typename MapType::tagType tagType;

  mappedVector(const MapType & _mapping): mapping(&_mapping), NumericVector<T>(_mapping.size()){}
  
  typedef mappedVector<T,MapType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,mappedVector<T,MapType> >::value, int>::type = 0>
  mappedVector(U&& expr){
    mapping = expr.common_properties();
    this->resize(mapping->size());
    for(int i=0;i<mapping->size();i++) this->operator()(i) = expr[i];
  }
  
  inline T & operator()(const tagType &tag){ return this->operator()(mapping->map(tag)); }
  inline const T & operator()(const tagType &tag) const{ return this->operator()(mapping->map(tag)); }
  inline T & operator()(const int idx){ return NumericVector<T>::operator()(idx); }
  inline const T & operator()(const int idx) const{ return NumericVector<T>::operator()(idx); }
  
  const MapType & getMapping() const{ return *mapping; }
};
template<typename T, typename MapType>
struct getElem<mappedVector<T,MapType> >{
  static inline auto elem(const mappedVector<T,MapType> &v, const int i)->decltype(v(i)){ return v(i); }
  // static inline auto elem(mappedVector<T,MapType> &v, const int i)->decltype(v(i)){ return v(i); }
  static inline MapType const* common_properties(const mappedVector<T,MapType> &v){ return &v.getMapping(); }
};

template<typename T, typename MapType>
inline void debug_print(const mappedVector<T,MapType> &v){ std::cout << &v.getMapping() << std::endl; std::cout.flush(); }

CPSFIT_END_NAMESPACE
#endif
