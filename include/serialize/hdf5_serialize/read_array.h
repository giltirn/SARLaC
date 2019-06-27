#ifndef _HDF5_READ_ARRAY_H___
#define _HDF5_READ_ARRAY_H___

//Functions to read array types
#include<config.h>

#ifdef HAVE_HDF5

#include<serialize/hdf5_serialize/hdf5_reader.h>
#include<serialize/hdf5_serialize/read_write_basic.h>

CPSFIT_START_NAMESPACE

//For (contiguous) arrays we can read in a compact format if native, otherwise we can read in an un-compact format. 
//Requires a policy indicating how to allocate and access (cf above)
template<typename ArrayPolicy, IF_NATIVE(typename ArrayPolicy::ElementType)>
inline static void readCompact(HDF5reader &reader, typename ArrayPolicy::ArrayType &v, const std::string &tag){
  reader.read<ArrayPolicy>(v,tag);
}

template<typename ArrayPolicy>
inline static void readUncompact(HDF5reader &reader, typename ArrayPolicy::ArrayType &value, const std::string &tag){
  reader.enter(tag); //enter a group
  unsigned long size;
  reader.read(size,"size");
  ArrayPolicy::resize(value, size);
  typename ArrayPolicy::ElementType *ptr = ArrayPolicy::pointer(value);

  for(unsigned long i=0;i<size;i++){
    std::ostringstream os;
    os << "elem_" << i;
    read(reader,ptr[i],os.str());
  }
  reader.leave();
}

//For all array types, those with native elements default to compact reads, and those with non-native to uncompact reads
template<typename ArrayPolicy, IF_NATIVE(typename ArrayPolicy::ElementType)>
inline static void read(HDF5reader &reader, typename ArrayPolicy::ArrayType &v, const std::string &tag){
  readCompact<ArrayPolicy>(reader,v,tag);
}
template<typename ArrayPolicy, IF_NOT_NATIVE(typename ArrayPolicy::ElementType)>
inline static void read(HDF5reader &reader, typename ArrayPolicy::ArrayType &v, const std::string &tag){
  readUncompact<ArrayPolicy>(reader,v,tag);
}


//C-arrays
template<typename T>
struct HDF5readerCArrayPolicy{
  typedef T* ArrayType;
  typedef T ElementType;
  inline static void resize(ArrayType &v, const hsize_t sz){ v = (T*)malloc(sz * sizeof(T)); }
  inline static ElementType* pointer(ArrayType &v){ return v; }
};
template<typename T>
inline static void readCompact(HDF5reader &reader, T* &v, const std::string &tag){
  readCompact<HDF5readerCArrayPolicy<T> >(reader, v, tag);
}
template<typename T>
inline static void readUncompact(HDF5reader &reader, T* &v, const std::string &tag){
  readUncompact<HDF5readerCArrayPolicy<T> >(reader, v, tag);
}
template<typename T>
inline static void read(HDF5reader &reader, T* &v, const std::string &tag){
  read<HDF5readerCArrayPolicy<T> >(reader, v, tag);
}

//Strings (only compact)
inline static void read(HDF5reader &reader, std::string &v, const std::string &tag){
  struct HDF5readerStringPolicy{
    typedef std::string ArrayType;
    typedef char ElementType;
    inline static void resize(ArrayType &v, const int sz){ v.resize(sz); }
    inline static ElementType* pointer(ArrayType &v){ return &v[0]; }
  };
  read<HDF5readerStringPolicy>(reader, v, tag);
}

//Complex (only compact)
template<typename T, IF_NATIVE(T)>
inline static void read(HDF5reader &reader, std::complex<T> &v, const std::string &tag){
  struct HDF5readerComplexPolicy{
    typedef std::complex<T> ArrayType;
    typedef T ElementType;
    inline static void resize(ArrayType &v, const int sz){ assert(sz == 2); }
    inline static ElementType* pointer(ArrayType &v){ return reinterpret_cast<T*>(&v); }
  };
  read<HDF5readerComplexPolicy>(reader, v, tag);
}

//std::vector
template<typename T>
struct HDF5readerVectorPolicy{
  typedef std::vector<T> ArrayType;
  typedef T ElementType;
  inline static void resize(ArrayType &v, const int sz){ v.resize(sz); }
  inline static ElementType* pointer(ArrayType &v){ return v.data(); }
};
template<typename T>
inline static void readCompact(HDF5reader &reader, std::vector<T> &v, const std::string &tag){
  readCompact<HDF5readerVectorPolicy<T> >(reader, v, tag);
}
template<typename T>
inline static void readUncompact(HDF5reader &reader, std::vector<T> &v, const std::string &tag){
  readUncompact<HDF5readerVectorPolicy<T> >(reader, v, tag);
}
//For distributions the defaults are specified elsewhere
template<typename T, IF_NOT_DISTRIBUTION_NATIVE(T)>
inline static void read(HDF5reader &reader, std::vector<T> &v, const std::string &tag){
  read<HDF5readerVectorPolicy<T> >(reader, v, tag);
}

//Overload compact reads for complex
template<typename T>
struct HDF5readerVectorComplexPolicy{
  typedef std::vector<std::complex<T> > ArrayType;
  typedef T ElementType;
  inline static void resize(ArrayType &v, const int sz){ assert(sz % 2 == 0); v.resize(sz/2); }
  inline static ElementType* pointer(ArrayType &v){ return reinterpret_cast<T*>(v.data()); }
};
template<typename T, IF_NATIVE(T)>
inline static void readCompact(HDF5reader &reader, std::vector<std::complex<T> > &v, const std::string &tag){
  readCompact<HDF5readerVectorComplexPolicy<T> >(reader, v, tag);
}
template<typename T, IF_NATIVE(T)>
inline static void read(HDF5reader &reader, std::vector<std::complex<T> > &v, const std::string &tag){
  //For backwards compatibility we test whether the array was written in uncompacted vector format, which was the old default
  //This can be done by checking if a group with name 'tag' exists (if written in compact format it is an attribute instead)
  if(reader.containsGroup(tag)) readUncompact<HDF5readerVectorPolicy<std::complex<T> > >(reader, v, tag);
  else readCompact<HDF5readerVectorComplexPolicy<T> >(reader, v, tag);
}


template<typename T, std::size_t Size>
struct HDF5readerArrayPolicy{
  typedef std::array<T,Size> ArrayType;
  typedef T ElementType;
  inline static void resize(ArrayType &v, const int sz){ assert(sz == Size); }
  inline static ElementType* pointer(ArrayType &v){ return v.data(); }
};
template<typename T, int size>
inline static void readCompact(HDF5reader &reader, std::array<T,size> &v, const std::string &tag){
  readCompact<HDF5readerArrayPolicy<T,size> >(reader, v, tag);
}
template<typename T, int size>
inline static void readUncompact(HDF5reader &reader, std::array<T,size> &v, const std::string &tag){
  readUncompact<HDF5readerArrayPolicy<T,size> >(reader, v, tag);
}
template<typename T, int size>
inline static void read(HDF5reader &reader, std::array<T,size> &v, const std::string &tag){
  read<HDF5readerArrayPolicy<T,size> >(reader, v, tag);
}

CPSFIT_END_NAMESPACE

#endif

#endif
