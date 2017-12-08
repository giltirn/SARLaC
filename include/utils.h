#ifndef _UTILS_H__
#define _UTILS_H__
#include<cstdarg>
#include<cassert>
#include<cstdio>
#include<memory>
#include<iostream>
#include<cxxabi.h>
#include<sstream>
#include<limits>
#include<omp.h>
#include<sys/sysinfo.h>

#include<template_wizardry.h>

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



template<typename T, typename std::enable_if<hasZeroMethod<T>::value, int>::type = 0>
inline void zeroit(T &v){
  v.zero();
}
template<typename T, typename std::enable_if<!hasZeroMethod<T>::value && (hasEqualsMethod<T,double>::value || hasEqualsMethod<T,float>::value), int>::type = 0>
inline void zeroit(T &v){
  v = 0.;
}
inline void zeroit(double &v){
  v = 0.;
}
inline void zeroit(float &v){
  v = 0.;
}

template<typename Operation>
auto threadedSum(const Operation &op)->typename std::decay<decltype(op(0))>::type{
  typedef typename std::decay<decltype(op(0))>::type T;
  const int N = op.size();
  const int nthread = omp_get_max_threads();
 
  T init_zero(op(0)); zeroit(init_zero);    //_threadedSumHelper<Operation,T>::getZero(op);
  std::vector<T> sum(nthread, init_zero);

#pragma omp parallel for
  for(int i=0;i<N;i++)
    sum[omp_get_thread_num()] = sum[omp_get_thread_num()] + op(i);

  for(int t=1;t<nthread;t++)
    sum[0] = sum[0] + sum[t];

  return sum[0];
}

template<typename T>
T threadedSum(const std::vector<T> &v){
  struct Op{
    const std::vector<T> &vv;
    inline const T & operator()(const int idx) const{ return vv[idx]; }
    inline size_t size() const{ return vv.size(); }
    Op(const std::vector<T> &_vv): vv(_vv){}
  };
  Op op(v);
  return threadedSum<Op>(op);
}

inline std::string stringize(const char* format, ...){
  int n; //not counting null character
  {
    char buf[1024];
    va_list argptr;
    va_start(argptr, format);    
    n = vsnprintf(buf, 1024, format, argptr);
    va_end(argptr);
    if(n < 1024) return std::string(buf);
  }
  
  char buf[n+1];
  va_list argptr;
  va_start(argptr, format);    
  int n2 = vsnprintf(buf, 1024, format, argptr);
  va_end(argptr);
  assert(n2 <= n);
  return std::string(buf);
}



//Print memory usage
inline void printMem(const std::string &reason = "", FILE* stream = stdout){
  if(reason != "") fprintf(stream, "printMem called with reason '%s': ", reason.c_str());
  else fprintf(stream, "printMem: ");
  
  struct sysinfo myinfo;
  sysinfo(&myinfo);  //cf http://man7.org/linux/man-pages/man2/sysinfo.2.html
  double total_mem = myinfo.mem_unit * myinfo.totalram;
  total_mem /= (1024.*1024.);
  double free_mem = myinfo.mem_unit * myinfo.freeram;
  free_mem /= (1024.*1024.);

  fprintf(stream,"Memory: total: %.2f MB, avail: %.2f MB, used %.2f MB\n",total_mem, free_mem, total_mem-free_mem);

#define REPORT_MEMINFO
#ifdef REPORT_MEMINFO
  {
    double total_mem2;
    double free_mem2;
  
    std::string token;
    std::ifstream file("/proc/meminfo");
    while(file >> token) {
        if(token == "MemTotal:") {
	  file >> total_mem2;
	}else if(token == "MemAvailable:"){
	  file >> free_mem2;
	}
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');        // ignore rest of the line
    }
    file.close();

    fprintf(stream,"Memory (/proc/meminfo): total: %.2f MB, avail: %.2f MB, used %.2f MB\n",total_mem2/1024., free_mem2/1024., (total_mem2-free_mem2)/1024.);
  }
#endif
  
  fflush(stream);
}

inline bool fileExists(const std::string &filename){
  std::ifstream infile(filename);
  return infile.good();
}

#endif
