#ifndef _HDF5_TYPE_MAP_H_
#define _HDF5_TYPE_MAP_H_

#include<config.h>

#ifdef HAVE_HDF5

#include<H5Cpp.h>

#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//Wrap the awkward HDF5 types in template specializations to define a compile-time mapping
template<typename T>
struct H5typeMap{
  enum {is_native = 0};
};

#define H5_TYPE_MAP(TYPE, TYPENAME)					\
  template<>								\
  struct H5typeMap<TYPE>{						\
    static inline const H5::DataType & type(void){ return H5::PredType::TYPENAME; } \
    enum {is_native = 1};						\
  }

H5_TYPE_MAP(bool, NATIVE_B8);
H5_TYPE_MAP(char, NATIVE_CHAR);
H5_TYPE_MAP(signed char, NATIVE_SCHAR);
H5_TYPE_MAP(unsigned char, NATIVE_UCHAR);
H5_TYPE_MAP(short, NATIVE_SHORT);
H5_TYPE_MAP(unsigned short, NATIVE_USHORT);
H5_TYPE_MAP(int, NATIVE_INT);
H5_TYPE_MAP(unsigned int, NATIVE_UINT);
H5_TYPE_MAP(long, NATIVE_LONG);
H5_TYPE_MAP(unsigned long, NATIVE_ULONG);
H5_TYPE_MAP(long long, NATIVE_LLONG);
H5_TYPE_MAP(unsigned long long, NATIVE_ULLONG);
H5_TYPE_MAP(float, NATIVE_FLOAT);
H5_TYPE_MAP(double, NATIVE_DOUBLE);
H5_TYPE_MAP(long double, NATIVE_LDOUBLE);

CPSFIT_END_NAMESPACE

#endif

#endif
