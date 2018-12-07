#ifndef _HDF5_READWRITE_BASIC_H___
#define _HDF5_READWRITE_BASIC_H___

//Functions to read and write basic types (POD types, vectors of POD types, strings)
#include<config.h>

#ifdef HAVE_HDF5

#include<serialize/hdf5_serialize/hdf5_writer.h>
#include<serialize/hdf5_serialize/hdf5_reader.h>

CPSFIT_START_NAMESPACE

//Native types and vectors thereof
template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline void write(HDF5writer &writer, const T &value, const std::string &tag){
  writer.write(value,tag);
}
template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline void read(HDF5reader &reader, T &value, const std::string &tag){
  reader.read(value,tag);
}
template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline static void write(HDF5writer &writer, const std::vector<T> &value, const std::string &tag){
  writer.write(value,tag);
}
template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline static void read(HDF5reader &reader, std::vector<T> &value, const std::string &tag){
  reader.read(value,tag);
}
template<typename T, std::size_t Size,typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline static void write(HDF5writer &writer, const std::array<T,Size> &value, const std::string &tag){
  writer.write(value,tag);
}
template<typename T, std::size_t Size,typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline static void read(HDF5reader &reader, std::array<T,Size> &value, const std::string &tag){
  reader.read(value,tag);
}
template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline static void write(HDF5writer &writer, const std::complex<T> &value, const std::string &tag){
  writer.write(value,tag);
}
template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
inline static void read(HDF5reader &reader, std::complex<T> &value, const std::string &tag){
  reader.read(value,tag);
}

//Strings
inline void write(HDF5writer &writer, const std::string &value, const std::string &tag){
  writer.write(value,tag);
}
inline void read(HDF5reader &reader, std::string &value, const std::string &tag){
  reader.read(value,tag);
}

CPSFIT_END_NAMESPACE

#endif

#endif
