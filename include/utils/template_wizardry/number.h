#ifndef _CPSFIT_TEMPLATE_WIZARDRY_NUMBER_H_
#define _CPSFIT_TEMPLATE_WIZARDRY_NUMBER_H_

//Metaprogramming constructs for obtaining information about numerical types

#include<config.h>
#include<utils/macros.h>
#include<complex>

CPSFIT_START_NAMESPACE

//Check if type is an std::complex
template<typename T>
struct is_std_complex{ enum {value = 0}; };

template<typename T>
struct is_std_complex<std::complex<T> >{ enum {value = 1}; };

//Check if a type is a real or imaginary scalar (double, float, int)
template<typename T>
struct is_scalar{ enum{value = 0}; };

template<>
struct is_scalar<double>{ enum{value = 1}; };

template<>
struct is_scalar<float>{ enum{value = 1}; };

template<>
struct is_scalar<int>{ enum{value = 1}; };

template<typename T>
struct is_scalar<std::complex<T> >{ enum{value = is_scalar<T>::value }; };

//Check if type is a real or imaginary number
template<typename T>
struct is_floating_point_or_complex{ enum { value = std::is_floating_point<T>::value || is_std_complex<T>::value }; };

CPSFIT_END_NAMESPACE
#endif
