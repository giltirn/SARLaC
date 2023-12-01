#ifndef _GENERIC_ET_LOOP_EVAL_H_
#define _GENERIC_ET_LOOP_EVAL_H_

#include<ET/generic_ET/value_storage.h>

SARLAC_START_NAMESPACE

template<typename T, typename U>
inline void loop_eval(T &obj, U &&expr){
    for(int i=0;i<obj.size();i++) ETeval<T>::elem(obj,i) = expr[i];
}

template<typename T, typename U>
inline void parallel_loop_eval(T &obj, U &&expr){
#pragma omp parallel for
    for(int i=0;i<obj.size();i++) ETeval<T>::elem(obj,i) = expr[i];
}



SARLAC_END_NAMESPACE

#endif
