#ifndef _CPSFIT_UTILS_TEXT_OUTPUT_
#define _CPSFIT_UTILS_TEXT_OUTPUT_

//Functions and types for ascii output
#include<vector>
#include<array>
#include<iostream>

#include<config.h>
#include<utils/macros.h>


CPSFIT_START_NAMESPACE

//Any derived type DerivedType of OstreamHook with a write(std::ostream &) method will automatically have an operator<<(std::ostream &os, const DerivedType &v)
class OstreamHook{
public:
  virtual void write(std::ostream &) const = 0;
};

inline std::ostream & operator<<(std::ostream &os, const OstreamHook &hk){
  hk.write(os);
  return os;
}

//Ostream output for vector
template<typename T>
std::ostream & operator<<(std::ostream &os, const std::vector<T> &s){
  os << '(';
  if(s.size() > 0){
    for(int i=0;i<s.size()-1;i++) os << s[i] << ", ";
    os << s.back();
  }
  os << ')';
  return os;
}

template<typename T, std::size_t S>
std::ostream & operator<<(std::ostream &os, const std::array<T,S> &s){
  os << '(';
  if(S > 0){
    for(int i=0;i<s.size()-1;i++) os << s[i] << ", ";
    os << s.back();
  }
  os << ')';
  return os;
}


//An ostream slot-in that just throws text away without printing - useful for suppressing output in threaded environments
struct nullstream: std::ostream {
  nullstream(): std::ios(0), std::ostream(0) {}
};
static nullstream null_stream;



CPSFIT_END_NAMESPACE
#endif
