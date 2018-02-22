#ifndef _HDF5_WRITER_H___
#define _HDF5_WRITER_H___

#include<config.h>

//A class for writing to HDF5 files

#ifdef HAVE_HDF5

#include<serialize/hdf5_serialize/type_map.h>
#include<utils/utils.h>

CPSFIT_START_NAMESPACE

class HDF5writer{
  H5::H5File file;
  std::vector<H5::Group> group;

  //Write array types
  template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
  void write(T const* v, const int sz, const std::string &name){
    using namespace H5;
    const DataType &type = H5typeMap<T>::type();
    hsize_t len(sz);
    DataSpace dataSpace(1, &len);
    std::size_t byte_size = sz * sizeof(T);    
    if(byte_size > 64*1024){ //write as a data set
      DataSet dataSet = group.back().createDataSet(name.c_str(),type, dataSpace);
      dataSet.write(v,type);
    }else{
      Attribute attribute = group.back().createAttribute(name.c_str(),type, dataSpace);
      attribute.write(type,v);
    }  
  }
public:
  HDF5writer(const std::string &filename): file(filename.c_str(), H5F_ACC_TRUNC ){
    H5::Exception::dontPrint();
    group.push_back(file.openGroup("/"));
  }

  //Write a single value as an attribute
  template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
  void write(const T &v, const std::string &name){
    using namespace H5;
    hsize_t attrDim = 1;
    DataSpace attrSpace(1, &attrDim);
    const DataType &type = H5typeMap<T>::type();
    Attribute attribute = group.back().createAttribute(name.c_str(), type, attrSpace);
    attribute.write(type, &v);    
  }

  //Write vector types
  template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
  inline void write(const std::vector<T> &v, const std::string &name){
    write(v.data(),v.size(),name);
  }
  //Write array types
  template<typename T, std::size_t Size, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
  inline void write(const std::array<T,Size> &v, const std::string &name){
    write(v.data(),Size,name);
  }
  //Write string
  void write(const std::string &v, const std::string &name){
    return write(v.data(),v.size(),name);
  }
  
  void enter(const std::string &nm){
    try{
      group.push_back(group.back().openGroup(nm.c_str()));
    }catch(H5::Exception& e){
      group.push_back(group.back().createGroup(nm.c_str()));
    }
  }
  void leave(){
    group.pop_back();
  }    
};



CPSFIT_END_NAMESPACE

#endif

#endif
