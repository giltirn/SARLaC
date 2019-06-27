#ifndef _HDF5_READWRITE_NONNATIVE_H___
#define _HDF5_READWRITE_NONNATIVE_H___

//Functions to read and write arrays of non-native types. The non-native type must have an implementation of read and write.
#include<config.h>

#ifdef HAVE_HDF5

#include<utils/template_wizardry.h>
#include<serialize/hdf5_serialize/read_write_basic.h>
#include<serialize/hdf5_serialize/read_array.h>
#include<serialize/hdf5_serialize/write_array.h>

CPSFIT_START_NAMESPACE

//Non-native vector<vector> collapse the two arrays into one to avoid metadata overheads
template<typename T, IF_NOT_NATIVE(T), IF_NOT_DISTRIBUTION_NATIVE(T)>
inline static void write(HDF5writer &writer, const std::vector<std::vector<T> > &value, const std::string &tag){
  writer.enter(tag); //enter a group
  unsigned long size1 = value.size();
  write(writer, size1,"size1");
  
  std::vector<unsigned long> size2(value.size());
  for(int i=0;i<value.size();i++) size2[i] = value[i].size();  
  writeCompact(writer, size2,"size2");
  
  for(int i=0;i<value.size();i++){
    for(int j=0;j<value[i].size();j++){    
      std::ostringstream os;
      os << "elem_" << i << "_" << j;
      write(writer,value[i][j],os.str());
    }
  }
  writer.leave();
}
template<typename T, IF_NOT_NATIVE(T), IF_NOT_DISTRIBUTION_NATIVE(T)>
inline static void read(HDF5reader &reader, std::vector<std::vector<T> > &value, const std::string &tag){
  reader.enter(tag); //enter a group
  unsigned long size1;
  read(reader, size1,"size1");
  
  std::vector<unsigned long> size2;
  readCompact(reader, size2,"size2");

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






CPSFIT_END_NAMESPACE

#endif

#endif
