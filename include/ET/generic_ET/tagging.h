#ifndef _GENERIC_ET_TAGGING_H_
#define _GENERIC_ET_TAGGING_H_

#include<utils/macros.h>
#include<utils/template_wizardry.h>

CPSFIT_START_NAMESPACE

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


CPSFIT_END_NAMESPACE

#endif
