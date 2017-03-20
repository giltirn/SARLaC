#ifndef _TEMPLATE_WIZARDRY_H_
#define _TEMPLATE_WIZARDRY_H_
//Now with added c++11!

#include<complex>
#include<type_traits>


template<typename T>
struct is_std_complex{ enum {value = 0}; };

template<typename T>
struct is_std_complex<std::complex<T> >{ enum {value = 1}; };

template<typename T>
struct is_std_vector{ enum {value = 0}; };

template<typename T>
struct is_std_vector<std::vector<T> >{ enum {value = 1}; };

template<typename T>
struct is_floating_point_or_complex{ enum { value = std::is_floating_point<T>::value || is_std_complex<T>::value }; };

template<typename T, typename U>
struct add_const_if{
  typedef T type;
};
template<typename T, typename U>
struct add_const_if<T, const U>{
  typedef const T type;
};

template<typename VectorType>
struct get_value_type{
  typedef typename std::remove_reference<
  decltype( ((VectorType*)(nullptr))->operator[](0) )
    >::type type;
};

template<typename VectorType, typename T>
struct value_type_equals{
  enum{ value = std::is_same< typename get_value_type<VectorType>::type, T >::value };
};

template<class T>
struct Void {
  typedef void type;
};

template<typename T, typename U = void>
struct hasSampleMethod{
  enum{ value = 0 };
};
template<typename T>
struct hasSampleMethod<T, typename Void<decltype( ((T*)(NULL))->sample(0) )>::type>{
  enum{ value = 1 };
};
template<typename T, typename U = void>
struct hasMeanMethod{
  enum{ value = 0 };
};
template<typename T>
struct hasMeanMethod<T, typename Void<decltype( ((T*)(NULL))->mean() )>::type>{
  enum{ value = 1 };
};

template<class T, class Fallback = void>
struct hasDataType{ enum{ value = 0 }; };

template<class T>
struct hasDataType<T, typename Void<typename T::DataType>::type >{ enum{ value = 1 }; };

#endif
