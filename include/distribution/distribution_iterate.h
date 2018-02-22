#ifndef _DISTRIBUTION_ITERATE_H
#define _DISTRIBUTION_ITERATE_H

#include<iostream>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>

CPSFIT_START_NAMESPACE

//This type of iterator is bound to a distribution instance
template<typename T_unconst, int is_const>
class _distributionIterator{
  typedef typename add_const_if_int<T_unconst,is_const>::type distributionType;
  typedef typename add_const_if_int<typename T_unconst::DataType,is_const>::type type;
  int s;
  int sz;
  distributionType* d;
public:
  _distributionIterator(): s(0), sz(0), d(NULL){}
  _distributionIterator(distributionType &dist): s(0), sz(dist.size()), d(&dist){}
  inline void operator++(){ ++s; }
  inline bool end() const{ return s>=sz; }
  inline type& operator*() const{ return d->sample(s); }
  inline void report() const{ std::cout << s << " (" << sz << ")" << std::endl; }
};

template<typename T> using distributionIterator = _distributionIterator<typename std::remove_const<T>::type, std::is_const<T>::value>;

//These structures allow you to generically iterate over the elements of many distributions (of the same size)
//Should be individually specialized for the different distribution types
template<typename distributionType>
struct iterate;


CPSFIT_END_NAMESPACE

#endif
