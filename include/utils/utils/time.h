#ifndef _SARLAC_UTILS_TIME_
#define _SARLAC_UTILS_TIME_

//Functions for timing
#include<chrono>

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

struct timer{
  std::chrono::time_point<std::chrono::high_resolution_clock> _start;
  std::chrono::time_point<std::chrono::high_resolution_clock> _stop;

  inline void start(){ _start = std::chrono::high_resolution_clock::now(); }
  inline void stop(){ _stop = std::chrono::high_resolution_clock::now(); }
  inline double elapsed(){ return std::chrono::duration<double, std::nano>(_stop-_start).count(); } //in nanoseconds
};

SARLAC_END_NAMESPACE

#endif
