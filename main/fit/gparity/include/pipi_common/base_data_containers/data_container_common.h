#pragma once

#include<map>
#include <boost/timer/timer.hpp>

#include<config.h>
#include<utils/macros.h>

#include<distribution.h>
#include<tensors.h>
#include<common.h>

SARLAC_START_NAMESPACE

enum SourceOrSink { Source, Sink };
inline void write(SARLaC::HDF5writer &writer, const SourceOrSink d, const std::string &tag){ write(writer,(int)d,tag); }
inline void read(SARLaC::HDF5reader &reader, SourceOrSink &d, const std::string &tag){
  int dd; read(reader,dd,tag); 
  d = (SourceOrSink)dd;
}

//Useful policies for setting output type to real or re/im
struct setValueRealPart{
  template<typename T>
  static inline void set(T &into, const double re, const double im){ into = re; }
};
struct setValueComplex{
  template<typename T>
  static inline void set(T &into, const double re, const double im){ into.real(re); into.imag(im); }
};

SARLAC_END_NAMESPACE
