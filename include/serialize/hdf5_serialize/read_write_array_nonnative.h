#ifndef _HDF5_READWRITE_NONNATIVE_H___
#define _HDF5_READWRITE_NONNATIVE_H___

//Functions to read and write arrays of non-native types. The non-native type must have an implementation of read and write.
#include<config.h>

#ifdef HAVE_HDF5

#include<utils/template_wizardry.h>
#include<serialize/hdf5_serialize/read_write_basic.h>

CPSFIT_START_NAMESPACE

template<typename T, int>
struct _isDistributionOfHDF5nativetype{ enum {value = 0}; };

template<typename T>
struct _isDistributionOfHDF5nativetype<T,1>{ enum {value = H5typeMap<typename T::DataType>::is_native}; }; 

template<typename T>
struct isDistributionOfHDF5nativetype{ enum {value = _isDistributionOfHDF5nativetype<T,hasSampleMethod<T>::value>::value }; };    

  
//Non-native vectors (but not distributions, they are specialized elsewhere)
template<typename T, typename std::enable_if<!H5typeMap<T>::is_native && !isDistributionOfHDF5nativetype<T>::value, int>::type = 0>
inline static void write(HDF5writer &writer, const std::vector<T> &value, const std::string &tag){
  writer.enter(tag); //enter a group
  unsigned long size = value.size();
  writer.write(size,"size");  
  for(int i=0;i<value.size();i++){
    std::ostringstream os;
    os << "elem_" << i;
    write(writer,value[i],os.str());
  }
  writer.leave();
}
template<typename T, typename std::enable_if<!H5typeMap<T>::is_native && !isDistributionOfHDF5nativetype<T>::value, int>::type = 0>
inline static void read(HDF5reader &reader, std::vector<T> &value, const std::string &tag){
  reader.enter(tag); //enter a group
  unsigned long size;
  reader.read(size,"size");
  value.resize(size);
  for(int i=0;i<value.size();i++){
    std::ostringstream os;
    os << "elem_" << i;
    read(reader,value[i],os.str());
  }
  reader.leave();
}

//Non-native arrays
template<typename T, std::size_t Size, typename std::enable_if<!H5typeMap<T>::is_native, int>::type = 0>
inline static void write(HDF5writer &writer, const std::array<T,Size> &value, const std::string &tag){
  writer.enter(tag); //enter a group
  unsigned long size = value.size();
  writer.write(size,"size");  
  for(int i=0;i<value.size();i++){
    std::ostringstream os;
    os << "elem_" << i;
    write(writer,value[i],os.str());
  }
  writer.leave();
}
template<typename T, std::size_t Size, typename std::enable_if<!H5typeMap<T>::is_native, int>::type = 0>
inline static void read(HDF5reader &reader, std::array<T,Size> &value, const std::string &tag){
  reader.enter(tag); //enter a group
  unsigned long size;
  reader.read(size,"size");
  assert(size == Size);
  for(int i=0;i<value.size();i++){
    std::ostringstream os;
    os << "elem_" << i;
    read(reader,value[i],os.str());
  }
  reader.leave();
}


//Non-native vector<vector>
template<typename T, typename std::enable_if<!H5typeMap<T>::is_native && !isDistributionOfHDF5nativetype<T>::value, int>::type = 0>
inline static void write(HDF5writer &writer, const std::vector<std::vector<T> > &value, const std::string &tag){
  writer.enter(tag); //enter a group
  unsigned long size1 = value.size();
  writer.write(size1,"size1");
  
  std::vector<unsigned long> size2(value.size());
  for(int i=0;i<value.size();i++) size2[i] = value[i].size();  
  writer.write(size2,"size2");
  
  for(int i=0;i<value.size();i++){
    for(int j=0;j<value[i].size();j++){    
      std::ostringstream os;
      os << "elem_" << i << "_" << j;
      write(writer,value[i][j],os.str());
    }
  }
  writer.leave();
}
template<typename T, typename std::enable_if<!H5typeMap<T>::is_native && !isDistributionOfHDF5nativetype<T>::value, int>::type = 0>
inline static void read(HDF5reader &reader, std::vector<std::vector<T> > &value, const std::string &tag){
  reader.enter(tag); //enter a group
  unsigned long size1;
  reader.read(size1,"size1");
  
  std::vector<unsigned long> size2;
  reader.read(size2,"size2");

  value.resize(size1);
  for(int i=0;i<value.size();i++){
    value[i].resize(size2[i]);
    for(int j=0;j<value[i].size();j++){    
      std::ostringstream os;
      os << "elem_" << i << "_" << j;
      read(reader,value[i][j],os.str());
    }
  }
  reader.leave();
}



//Non-native complex
template<typename T, typename std::enable_if<!H5typeMap<T>::is_native, int>::type = 0>
inline static void write(HDF5writer &writer, const std::complex<T> &value, const std::string &tag){
  writer.enter(tag); //enter a group
  write(writer, value.real(), "real");
  write(writer, value.imag(), "imag");
  writer.leave();
}
template<typename T, typename std::enable_if<!H5typeMap<T>::is_native, int>::type = 0>
inline static void read(HDF5reader &reader, std::complex<T> &value, const std::string &tag){
  reader.enter(tag); //enter a group
  read(reader,reinterpret_cast<T(&)[2]>(value)[0],"real");
  read(reader,reinterpret_cast<T(&)[2]>(value)[1],"imag");
  reader.leave();
}


CPSFIT_END_NAMESPACE

#endif

#endif
