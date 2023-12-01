#ifndef _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_IO_FORMAT_H_
#define _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_IO_FORMAT_H_

//Classes that control the implementation of the standard format for different distribution types

#include<config.h>
#include<utils/macros.h>
#include<distribution/distribution_hdf5io_basic.h>
#include<distribution/distribution_hdf5io_conventional/utils.h>

SARLAC_START_NAMESPACE

//When the int is 1, typedef the underlying data type
template<typename T, int>
struct getDataType{};

template<typename T>
struct getDataType<T,0>{ typedef void type; };

template<typename T>
struct getDataType<T,1>{ typedef typename T::DataType type; };

#define IO_ENABLE_IF_POD(D) typename std::enable_if< hasSampleMethod<D>::value && std::is_arithmetic<typename getDataType<D,hasDataType<D>::value>::type>::value, int>::type = 0
#define IO_ENABLE_IF_NOT_POD(D) typename std::enable_if<hasSampleMethod<D>::value && !std::is_arithmetic<typename getDataType<D,hasDataType<D>::value>::type>::value, int>::type = 0

//Define how to write vectors and vector-vectors of distributions of *POD types*. By default we use the basic format without flattening (uncompressed)
template<typename T>
struct standardIOformat{
  inline static void write(HDF5writer &writer, const std::vector<T> &value, const std::string &tag){ SARLaC::write(writer,value,tag,false); } //no flattening
  inline static void read(HDF5reader &reader, std::vector<T> &value, const std::string &tag){ SARLaC::read(reader,value,tag,false); }
  inline static void write(HDF5writer &writer, const std::vector<std::vector<T> > &value, const std::string &tag){ SARLaC::write(writer,value,tag,false); } //no flattening
  inline static void read(HDF5reader &reader, std::vector<std::vector<T> > &value, const std::string &tag){ SARLaC::read(reader,value,tag,false); }  
};

SARLAC_END_NAMESPACE
#endif
