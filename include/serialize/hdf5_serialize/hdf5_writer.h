#ifndef _HDF5_WRITER_H___
#define _HDF5_WRITER_H___

#include<config.h>

//A class for writing to HDF5 files

#ifdef HAVE_HDF5

#include<complex>
#include<serialize/hdf5_serialize/type_map.h>
#include<utils/utils.h>

CPSFIT_START_NAMESPACE

class HDF5writer{
  H5::H5File file;
  std::vector<H5::Group> group;

public:
  HDF5writer(const std::string &filename): file(filename.c_str(), H5F_ACC_TRUNC ){
    H5::Exception::dontPrint();
    group.push_back(file.openGroup("/"));
  }

  //Write a single value as an attribute
  template<typename T, IF_NATIVE(T)>
  void write(const T &v, const std::string &name){
    using namespace H5;
    hsize_t attrDim = 1;
    DataSpace attrSpace(1, &attrDim);
    const DataType &type = H5typeMap<T>::type();
    Attribute attribute = group.back().createAttribute(name.c_str(), type, attrSpace);
    attribute.write(type, &v);    
  }

  //Write a generic array. Array type must be able to provide a pointer to a contiguous memory region.
  //Requires a policy class specifying:
  //ArrayType, ElementType (must be a native type),  static int size(const ArrayType &),  static ElementType const* pointer(const ArrayType &v)
  template<typename ArrayPolicy, IF_NATIVE(typename ArrayPolicy::ElementType)>
  void write(const typename ArrayPolicy::ArrayType &v, const std::string &name){
    typedef typename ArrayPolicy::ElementType T;
    using namespace H5;
    const DataType &type = H5typeMap<T>::type();
    hsize_t len = ArrayPolicy::size(v);
    DataSpace dataSpace(1, &len);
    std::size_t byte_size = len * sizeof(T);    
    T const * ptr = ArrayPolicy::pointer(v);
    if(byte_size > 64*1024){ //write as a data set
      DataSet dataSet = group.back().createDataSet(name.c_str(),type, dataSpace);
      dataSet.write(ptr,type);
    }else{
      Attribute attribute = group.back().createAttribute(name.c_str(),type, dataSpace);
      attribute.write(type,ptr);
    }  
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
