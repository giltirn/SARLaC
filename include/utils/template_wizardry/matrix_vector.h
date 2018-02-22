#ifndef _CPSFIT_TEMPLATE_WIZARDRY_MATRIX_VECTOR_H_
#define _CPSFIT_TEMPLATE_WIZARDRY_MATRIX_VECTOR_H_

//Metaprogramming constructs for obtaining information about arrays, vectors and matrices
#include<cstdlib>
#include<type_traits>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry/types.h>

CPSFIT_START_NAMESPACE

//Check if type is an std::vector
template<typename T>
struct is_std_vector{ enum {value = 0}; };

template<typename T>
struct is_std_vector<std::vector<T> >{ enum {value = 1}; };

//For array classes with operator[], determine the data type
template<typename VectorType>
struct get_value_type{
  typedef typename std::remove_reference<
  decltype( ((VectorType*)(nullptr))->operator[](0) )
    >::type type;
};

//Check if the data type of an array VectorType is equal to T
template<typename VectorType, typename T>
struct value_type_equals{
  enum{ value = std::is_same< typename get_value_type<VectorType>::type, T >::value };
};

//Check if a class T is matrix-like, i.e. has an operator()(int, int)
template<class T, class Fallback = void>
struct isMatrixType{ enum{ value = 0 }; };

template<class T>
struct isMatrixType<T, typename Void<decltype( ((T*)(NULL))->operator()(0,0) )  >::type >{ enum{ value = 1 }; };

//Get the element type of a matrix-like class, i.e. has an operator()(int, int)
template<typename MatrixType>
struct _get_elem_type{
  typedef decltype( ( (MatrixType*)(nullptr) )->operator()(0,0) ) RefType;
  typedef typename std::remove_reference<RefType>::type BaseType;
  typedef typename std::remove_const<BaseType>::type type;
};

//Get the element type of a vector-like class, i.e. has an operator()(int)
template<typename VectorType>
struct _get_vector_elem_type{
  typedef decltype( ( (VectorType*)(nullptr) )->operator()(0) ) RefType;
  typedef typename std::remove_reference<RefType>::type BaseType;
  typedef typename std::remove_const<BaseType>::type type;
};

//An enable-if directive for asserting a matrix-like class is a floating-point type
#define ENABLE_IF_ELEM_TYPE_FLOATINGPT(MatrixType)   typename std::enable_if< \
	   std::is_floating_point<typename _get_elem_type<MatrixType>::type>::value \
	   , int>::type = 0



//Check if a class T has a method resize(int)
template<typename T, typename U = void>
struct hasResizeMethod{
  enum{ value = 0 };
};
template<typename T>
struct hasResizeMethod<T, typename Void<decltype( ((T*)(NULL))->resize(0) )>::type>{
  enum{ value = 1 };
};

//Check if a class T has a method operator()(int) const
template<typename T, typename U = void>
struct hasParenthesesConstAccessor{
  enum{ value = 0 };
};
template<typename T>
struct hasParenthesesConstAccessor<T, typename Void<decltype( ((T const*)(NULL))->operator()(0) )>::type>{
  enum{ value = 1 };
};

//Check if a class T has a method operator[](int) const
template<typename T, typename U = void>
struct hasSquareBracketsConstAccessor{
  enum{ value = 0 };
};
template<typename T>
struct hasSquareBracketsConstAccessor<T, typename Void<decltype( ((T const*)(NULL))->operator[](0) )>::type>{
  enum{ value = 1 };
};

//Check if a class T has a method size()
template<typename T, typename U = void>
struct hasSizeMethod{
  enum{ value = 0 };
};
template<typename T>
struct hasSizeMethod<T, typename Void<decltype( ((T*)(NULL))->size() )>::type>{
  enum{ value = 1 };
};


CPSFIT_END_NAMESPACE
#endif
