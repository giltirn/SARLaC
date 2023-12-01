#ifndef _SARLAC_UTILS_ERROR_H_
#define _SARLAC_UTILS_ERROR_H_

#include<iostream>

#include<config.h>
#include<utils/macros.h>


SARLAC_START_NAMESPACE

//Flush the error message on the ostream before exiting with code
inline void error_exit(std::ostream &msg, const int code = -1){
  msg.flush();
  exit(code);
}

SARLAC_END_NAMESPACE
#endif
