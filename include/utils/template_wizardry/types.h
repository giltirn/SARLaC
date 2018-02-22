#ifndef _CPSFIT_TEMPLATE_WIZARDRY_TYPES_H_
#define _CPSFIT_TEMPLATE_WIZARDRY_TYPES_H_

//Some useful generic types

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//A null class
struct empty_t{};

template<class T>
struct Void {
  typedef void type;
};

CPSFIT_END_NAMESPACE
#endif
