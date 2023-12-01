#ifndef _UTILS_ENDIAN_H__
#define _UTILS_ENDIAN_H__
#include<boost/endian/conversion.hpp> 

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

inline bool isBigEndian(void){
    union {
        uint32_t i;
        unsigned char c[4];
    } bint = {0x01020304};
    return bint.c[0] == 1; 
}

template<int> struct _intMap{};
template<> struct _intMap<8>{ typedef uint64_t type; };
template<> struct _intMap<4>{ typedef uint32_t type; };

template<typename T>
inline T reverseEndianness(const T in){
  typedef typename _intMap<sizeof(T)>::type int_type;
  union{
    T i;
    int_type o;
  } bint;
  bint.i = in;
  int_type asint = boost::endian::endian_reverse(bint.o);
  union{
    int_type i;
    T o;
  } rint;
  rint.i = asint;
  return rint.o;
}

SARLAC_END_NAMESPACE

#endif
