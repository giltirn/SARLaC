#ifndef _GENERIC_ET_LOOP_EVAL_H_
#define _GENERIC_ET_LOOP_EVAL_H_

#include<ET/generic_ET/value_storage.h>

CPSFIT_START_NAMESPACE

//PARALLELIZE_DISTRIBUTION_ET disabled by default because it makes performance worse in most cases (better with very larger, O(10000) configs)
#ifdef PARALLELIZE_DISTRIBUTION_ET
template<typename T, typename U>
inline void loop_eval(T &obj, U &&expr){
#pragma omp parallel for
    for(int i=0;i<obj.size();i++) ETeval<T>::elem(obj,i) = expr[i];
}
#else
template<typename T, typename U>
inline void loop_eval(T &obj, U &&expr){
    for(int i=0;i<obj.size();i++) ETeval<T>::elem(obj,i) = expr[i];
}
#endif

CPSFIT_END_NAMESPACE

#endif
