#ifndef _COMPLEXIFY_TYPE_H_
#define _COMPLEXIFY_TYPE_H_

#include<config.h>
#include<utils/macros.h>
#include<complex>

#include "distribution.h"
#include "number.h"
#include "type_classify.h"

CPSFIT_START_NAMESPACE

//Determine the (std::)complex type corresponding to a scalar type that can be POD or a distribution

//Implementation stuff
namespace _get_complex_type_n{
  struct scalar_mark{};
  struct distribution_mark{};
  
  template<typename T, typename Mark>
  struct getComplex;
  
  template<typename T>
  struct getComplex<T, scalar_mark>{
    typedef std::complex<T> type;
  };

  template<typename T>
  struct getComplex<T, distribution_mark>{
    typedef typename getBaseType<T>::type baseType;
    typedef typename T::template rebase<std::complex<baseType> > type;
  };

  template<typename T, typename Mark>
  struct getReal;
  
  template<typename T>
  struct getReal<T, scalar_mark>{
    typedef typename T::value_type type;
  };

  template<typename T>
  struct getReal<T, distribution_mark>{
    typedef typename getBaseType<T>::type baseType;
    typedef typename T::template rebase<typename baseType::value_type> type;
  };

  template<typename T>
  using classify = typename TestElem< is_scalar<T>::value, scalar_mark,
				      TestElem< hasSampleMethod<T>::value, distribution_mark,
						LastElem
						>
				      >::type;
};
  

//eg Complexify<double> -> std::complex<double>
//   Complexify<jackknifeDistribution<double> > -> jackknifeDistribution<std::complex<double> >
template<typename T>
using Complexify = 
  typename _get_complex_type_n::getComplex<T, _get_complex_type_n::classify<T> >::type;

//Get the real type
//Assuming T = std::complex<U> or other complex type with value_type defined, or a distribution of such a type
template<typename T>
using Realify = 
  typename _get_complex_type_n::getReal<T, _get_complex_type_n::classify<T> >::type;


CPSFIT_END_NAMESPACE

#endif
