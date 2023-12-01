#ifndef _HDF5_READWRITE_BASIC_H___
#define _HDF5_READWRITE_BASIC_H___

//Functions to read and write HDF5 native types
#include<config.h>

#ifdef HAVE_HDF5

#include<serialize/hdf5_serialize/hdf5_writer.h>
#include<serialize/hdf5_serialize/hdf5_reader.h>

SARLAC_START_NAMESPACE

//Read native type (non-native should be user-specified unless defined elsewhere in library
template<typename T, IF_NATIVE(T)>
inline void read(HDF5reader &reader, T &value, const std::string &tag){
  reader.read(value,tag);
}

//Write native types (non-native types should be user-defined unless defined elsewhere in library)
template<typename T, IF_NATIVE(T)>
inline void write(HDF5writer &writer, const T &value, const std::string &tag){
  writer.write(value,tag);
}

SARLAC_END_NAMESPACE

#endif

#endif
