#ifndef _HDF5_WRITE_ARRAY_H___
#define _HDF5_WRITE_ARRAY_H___

//Functions to write array types
#include<config.h>

#ifdef HAVE_HDF5

#include<serialize/hdf5_serialize/hdf5_writer.h>
#include<serialize/hdf5_serialize/read_write_basic.h>

CPSFIT_START_NAMESPACE

//For arrays we can write in a compact format if native, otherwise we can write in an un-compact format
template<typename ArrayPolicy, IF_NATIVE(typename ArrayPolicy::ElementType)>
inline void writeCompact(HDF5writer &writer, const typename ArrayPolicy::ArrayType &value, const std::string &tag){
  writer.write<ArrayPolicy>(value,tag);
}

//Write an array of anything in a non-compact format
template<typename ArrayPolicy>
inline static void writeUncompact(HDF5writer &writer, const typename ArrayPolicy::ArrayType &value, const std::string &tag){
  writer.enter(tag); //enter a group
  unsigned long size = ArrayPolicy::size(value);
  typename ArrayPolicy::ElementType const* ptr = ArrayPolicy::pointer(value);
  writer.write(size,"size");  
  for(unsigned long i=0;i<size;i++){
    std::ostringstream os;
    os << "elem_" << i;
    write(writer,ptr[i],os.str());
  }
  writer.leave();
}

//Write of native defaults to compact and non-native to uncompact
template<typename ArrayPolicy, IF_NATIVE(typename ArrayPolicy::ElementType)>
inline void write(HDF5writer &writer, const typename ArrayPolicy::ArrayType &value, const std::string &tag){
  CPSfit::writeCompact<ArrayPolicy>(writer,value, tag);
}
template<typename ArrayPolicy, IF_NOT_NATIVE(typename ArrayPolicy::ElementType)>
inline void write(HDF5writer &writer, const typename ArrayPolicy::ArrayType &value, const std::string &tag){
  CPSfit::writeUncompact<ArrayPolicy>(writer,value, tag);
}

//C-arrays
template<typename T>
struct HDF5writerCArrayPolicy{
  struct ArrayCon{
    T const* ptr;
    hsize_t size;
    ArrayCon(T const* ptr, hsize_t size): ptr(ptr), size(size){}
  };
  typedef ArrayCon ArrayType;
  typedef T ElementType;

  inline static hsize_t size(const ArrayType &v){ return v.size; }
  inline static ElementType const* pointer(const ArrayType &v){ return v.ptr; }
};

template<typename T>
inline void writeCompact(HDF5writer &writer, T const* value, const size_t size, const std::string &tag){
  typename HDF5writerCArrayPolicy<T>::ArrayCon c(value,size);
  CPSfit::writeCompact<HDF5writerCArrayPolicy<T> >(writer, c, tag);
}
template<typename T>
inline void writeUncompact(HDF5writer &writer, T const* value, const size_t size, const std::string &tag){
  typename HDF5writerCArrayPolicy<T>::ArrayCon c(value,size);
  CPSfit::writeUncompact<HDF5writerCArrayPolicy<T> >(writer, c, tag);
}
template<typename T>
inline void write(HDF5writer &writer, T const* value, const size_t size, const std::string &tag){
  typename HDF5writerCArrayPolicy<T>::ArrayCon c(value,size);
  CPSfit::write<HDF5writerCArrayPolicy<T> >(writer, c, tag);
}


//Strings (only compact)
inline void write(HDF5writer &writer, const std::string &value, const std::string &tag){
  CPSfit::writeCompact(writer, value.data(), value.size(), tag);
}

//Complex (only compact)
template<typename T, IF_NATIVE(T)>
inline static void write(HDF5writer &writer, const std::complex<T> &value, const std::string &tag){
  CPSfit::writeCompact(writer, reinterpret_cast<T const*>(&value), 2, tag);
}

//std::vector
template<typename T>
inline static void writeCompact(HDF5writer &writer, const std::vector<T> &value, const std::string &tag){
  CPSfit::writeCompact(writer, value.data(), value.size(),tag);
}
template<typename T>
inline static void writeUncompact(HDF5writer &writer, const std::vector<T> &value, const std::string &tag){
  CPSfit::writeUncompact(writer, value.data(), value.size(),tag);
}
//For distributions these are defined elsewhere
template<typename T, IF_NOT_DISTRIBUTION_NATIVE(T)>
inline static void write(HDF5writer &writer, const std::vector<T> &value, const std::string &tag){
  CPSfit::write(writer, value.data(), value.size(),tag);
}


//Overload compact writes for complex
template<typename T, IF_NATIVE(T)>
inline static void writeCompact(HDF5writer &writer, const std::vector<std::complex<T> > &value, const std::string &tag){
  CPSfit::writeCompact(writer, reinterpret_cast<T const*>(value.data()), 2*value.size(),tag);
}
template<typename T, IF_NATIVE(T)>
inline static void write(HDF5writer &writer, const std::vector<std::complex<T> > &value, const std::string &tag){
  CPSfit::writeCompact(writer, reinterpret_cast<T const*>(value.data()), 2*value.size(),tag);
}



//std::array
template<typename T, size_t size>
inline static void writeCompact(HDF5writer &writer, const std::array<T,size> &value, const std::string &tag){
  CPSfit::writeCompact(writer, value.data(), size,tag);
}
template<typename T, size_t size>
inline static void writeUncompact(HDF5writer &writer, const std::array<T,size> &value, const std::string &tag){
  CPSfit::writeUncompact(writer, value.data(), value.size(),tag);
}
template<typename T, size_t size>
inline static void write(HDF5writer &writer, const std::array<T,size> &value, const std::string &tag){
  CPSfit::write(writer, value.data(), size,tag);
}



//Overload compact writes for complex
template<typename T, size_t size, IF_NATIVE(T)>
inline static void writeCompact(HDF5writer &writer, const std::array<std::complex<T>,size > &value, const std::string &tag){
  CPSfit::writeCompact(writer, reinterpret_cast<T const*>(value.data()), 2*size,tag);
}
template<typename T, size_t size, IF_NATIVE(T)>
inline static void write(HDF5writer &writer, const std::array<std::complex<T>, size> &value, const std::string &tag){
  CPSfit::writeCompact(writer, reinterpret_cast<T const*>(value.data()), 2*size,tag);
}


CPSFIT_END_NAMESPACE

#endif

#endif
