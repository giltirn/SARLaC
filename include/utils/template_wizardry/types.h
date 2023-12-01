#ifndef _SARLAC_TEMPLATE_WIZARDRY_TYPES_H_
#define _SARLAC_TEMPLATE_WIZARDRY_TYPES_H_

//Some useful generic types

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

//A null class
struct empty_t{};

template<class T>
struct Void {
  typedef void type;
};

SARLAC_END_NAMESPACE
#endif
