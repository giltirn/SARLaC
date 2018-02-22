#ifndef _CPSFIT_UTILS_TYPE_INFO_H_
#define _CPSFIT_UTILS_TYPE_INFO_H_
#include<cxxabi.h>
#include<string>
#include<memory>

#include<config.h>
#include<utils/macros.h>


CPSFIT_START_NAMESPACE

//Remove the name mangling
inline std::string demangle( const char* mangled_name ) {
  std::size_t len = 0 ;
  int status = 0 ;
  std::unique_ptr< char, decltype(&std::free) > ptr(__cxxabiv1::__cxa_demangle( mangled_name, nullptr, &len, &status ), &std::free ) ;
  return ptr.get();
}

//Print to string the type
template<typename T>
inline std::string printType(){ return demangle(typeid(T).name()); }


CPSFIT_END_NAMESPACE
#endif
