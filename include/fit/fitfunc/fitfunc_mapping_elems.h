#ifndef _CPSFIT_FITFUNC_MAPPING_ELEMS_H_
#define _CPSFIT_FITFUNC_MAPPING_ELEMS_H_

#include<cstdlib>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//Base type of elements
template<typename DataStructTo, typename DataStructFrom>
struct ParameterMapBase{
  virtual void copy(DataStructTo &to, const DataStructFrom &from) const = 0;
  virtual ParameterMapBase<DataStructTo,DataStructFrom> *clone() const = 0;
};

//Wrapper to safely contain, copy, delete elements stored by pointer to base
template<typename DataStructTo, typename DataStructFrom>
struct ParameterMapBaseWrap{
  ParameterMapBase<DataStructTo,DataStructFrom> *e;
public:
  inline void copy(DataStructTo &to, const DataStructFrom &from) const{ return e->copy(to,from); }
  ParameterMapBaseWrap(): e(NULL){}  
  ParameterMapBaseWrap(const ParameterMapBase<DataStructTo,DataStructFrom> &v): e(v.clone()){}
  ParameterMapBaseWrap(const ParameterMapBaseWrap<DataStructTo,DataStructFrom> &v): e(v.e->clone()){}
  
  ParameterMapBaseWrap &operator=(const ParameterMapBase<DataStructTo,DataStructFrom> &v){
    if(e!=NULL) delete e;
    e = v.clone();
  }
  
  ~ParameterMapBaseWrap(){ if(e!=NULL) delete e; }
};

//General mappings
template<typename T, typename DataStructTo, typename DataStructFrom>
struct ParameterMap: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toptr;
  T DataStructFrom::* fromptr;
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to.*toptr = from.*fromptr;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ParameterMap<T,DataStructTo, DataStructFrom>(*this); }
};
template<typename T, typename DataStructTo, typename DataStructFrom>
struct ParameterFreeze: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toptr;
  T value;
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to.*toptr = value;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ParameterFreeze<T,DataStructTo, DataStructFrom>(*this); }
};   
template<typename T, typename DataStructTo, typename DataStructFrom>
struct MemberArrayElementMap: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toarrayptr;
  T DataStructFrom::* fromarrayptr;
  int to_elem;
  int from_elem;
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    (to.*toarrayptr)[to_elem] = (from.*fromarrayptr)[from_elem];
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new MemberArrayElementMap<T,DataStructTo, DataStructFrom>(*this); }
};
template<typename T, typename DataStructTo, typename DataStructFrom>
struct MemberArrayElementFreeze: public ParameterMapBase<DataStructTo, DataStructFrom>{
  T DataStructTo::* toarrayptr;
  int to_elem;
  T value;
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    (to.*toarrayptr)[to_elem] = value;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new MemberArrayElementFreeze<T,DataStructTo, DataStructFrom>(*this); }
};
template<typename DataStructTo, typename DataStructFrom>
struct ArrayElementMap: public ParameterMapBase<DataStructTo, DataStructFrom>{
  int to_elem;
  int from_elem;

  ArrayElementMap(const int _to_elem, const int _from_elem): to_elem(_to_elem), from_elem(_from_elem){}
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to[to_elem] = from[from_elem];
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ArrayElementMap<DataStructTo, DataStructFrom>(*this); }
};
template<typename T, typename DataStructTo, typename DataStructFrom>
struct ArrayElementFreeze: public ParameterMapBase<DataStructTo, DataStructFrom>{
  int to_elem;
  T value;

  ArrayElementFreeze(const int _to_elem, const T &to): to_elem(_to_elem), value(to){}
  
  void copy(DataStructTo &to, const DataStructFrom &from) const{
    to[to_elem] = value;
  }
  ParameterMapBase<DataStructTo, DataStructFrom>* clone() const{ return new ArrayElementFreeze<T,DataStructTo, DataStructFrom>(*this); }
};


CPSFIT_END_NAMESPACE
#endif
