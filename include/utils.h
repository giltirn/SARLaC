#ifndef _UTILS_H__
#define _UTILS_H__

#include <memory>
#include <cxxabi.h>
#include<sstream>
#include<omp.h>

class OstreamHook{
public:
  virtual void write(std::ostream &) const = 0;
};

inline std::ostream & operator<<(std::ostream &os, const OstreamHook &hk){
  hk.write(os);
  return os;
}

template<typename T>
inline T & real_ref(std::complex<T> &z){
  return reinterpret_cast<T(&)[2]>(z)[0];
}
template<typename T>
inline T & imag_ref(std::complex<T> &z){
  return reinterpret_cast<T(&)[2]>(z)[1];
}

//Substitute substring '%d' with configuration idx
inline std::string subsIdx(const std::string fmt, const int idx){
  std::string::size_type off = fmt.find("%d");
  if(off == std::string::npos){
    std::cout << "Could not find substring \"%d\" in format string " << fmt << std::endl;
    std::cout.flush();
    exit(-1);
  }
  std::ostringstream os; os << idx;
  std::string out(fmt);
  out.replace(off,2,os.str());
  return out;
}

std::string demangle( const char* mangled_name ) {

  std::size_t len = 0 ;
  int status = 0 ;
  std::unique_ptr< char, decltype(&std::free) > ptr(
						    __cxxabiv1::__cxa_demangle( mangled_name, nullptr, &len, &status ), &std::free ) ;
  return ptr.get() ;
}

template<typename T>
inline std::string printType(){ return demangle(typeid(T).name()); }

inline void error_exit(std::ostream &msg, const int code = -1){
  msg.flush();
  exit(code);
}

template<typename T>
inline T strToAny(const std::string &str){
  T out;
  std::stringstream os(str); os >> out;
  return out;
}
template<typename T>
inline std::string anyToStr(const T &p){
  std::ostringstream os; os << p; return os.str();
}

#endif
