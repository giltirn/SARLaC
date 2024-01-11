#ifndef _SARLAC_UTILS_ERROR_H_
#define _SARLAC_UTILS_ERROR_H_

#include<iostream>

#include<config.h>
#include<utils/macros.h>
#include<gsl/gsl_errno.h>

SARLAC_START_NAMESPACE

//Flush the error message on the ostream before exiting with code
inline void error_exit(std::ostream &msg, const int code = -1){
  msg.flush();
  throw std::runtime_error(std::to_string(code));
}

inline void gsl_error_handler_warn(const char * reason, const char * file, int line, int gsl_errno){
  std::cout << "WARNING: GSL reported an error \"" << reason << "\" in file " << file << " line " << line << " with error code " << gsl_errno << std::endl;
}
inline void gsl_error_handler_throw(const char * reason, const char * file, int line, int gsl_errno){
  std::stringstream ss;
  ss << "GSL reported an error \"" << reason << "\" in file " << file << " line " << line << " with error code " << gsl_errno;
  throw std::runtime_error(ss.str());
}


inline void setGSLerrorHandlerWarn(){
  gsl_set_error_handler(&gsl_error_handler_warn);
}
inline void setGSLerrorHandlerThrow(){
  gsl_set_error_handler(&gsl_error_handler_throw);
}
inline void setGSLerrorHandlerAbort(){ //this is the default
  gsl_set_error_handler (NULL);
}



SARLAC_END_NAMESPACE
#endif
