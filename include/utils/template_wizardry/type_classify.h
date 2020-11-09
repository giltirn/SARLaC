#ifndef _TYPE_CLASSIFY_H_
#define _TYPE_CLASSIFY_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//A method of providing a list of conditions and associated types for classification

/*Usage: 
  template<typename T>
  using classify = typename TestElem< condition1<T>, truetype1,
                             TestElem< condition2<T>, truetype2,
	                               LastElem
                                     >
		           >::type; 
  Then if condition1<T> == true:
  typename classify<T>::type == truetype1
  
  if condition1<T> == false && condition2<T> == true
  typename classify<T>::type == truetype2
 */

struct no_mark{};

template<bool Condition, typename IfTrue, typename NextTest>
struct TestElem{};

template<typename IfTrue, typename NextTest>
struct TestElem<true,IfTrue,NextTest>{
  typedef IfTrue type;
};
template<typename IfTrue, typename NextTest>
struct TestElem<false,IfTrue,NextTest>{
  typedef typename NextTest::type type;
};

struct LastElem{
  typedef no_mark type;
};

/*
An if-else version of the above
Usage:
  template<typename T>
  using classify = typename TypeIfElse< condition<T>, truetype, falsetype >::type
*/

template<bool Condition, typename IfTrue, typename IfFalse>
struct TypeIfElse{};

template<typename IfTrue, typename IfFalse>
struct TypeIfElse<true,IfTrue,IfFalse>{
  typedef IfTrue type;
};

template<typename IfTrue, typename IfFalse>
struct TypeIfElse<false,IfTrue,IfFalse>{
  typedef IfFalse type;
};



CPSFIT_END_NAMESPACE

#endif
