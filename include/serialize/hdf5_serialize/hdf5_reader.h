#ifndef _HDF5_READER_H___
#define _HDF5_READER_H___

#include<config.h>

//A class for reading from HDF5 files

#ifdef HAVE_HDF5
#include<complex>
#include<serialize/hdf5_serialize/type_map.h>
#include<utils/utils.h>
CPSFIT_START_NAMESPACE

 
class HDF5reader{
  H5::H5File file;
  std::vector<H5::Group> group;

  static inline std::string filenameCheck(const std::string &filename){
    if(!fileExists(filename)) error_exit(std::cout << "HDF5reader constructor: File " << filename << " does not exist\n");
    return filename;
  }

public:
  HDF5reader(const std::string &filename): file(filenameCheck(filename).c_str(), H5F_ACC_RDONLY){
    H5::Exception::dontPrint();
    group.push_back(file.openGroup("/"));
  }

  //Read a single value as an attribute
  template<typename T, IF_NATIVE(T)>
  void read(T &v, const std::string &name){
    using namespace H5;
    Attribute attribute = group.back().openAttribute(name.c_str());
    DataSpace space = attribute.getSpace();
    hsize_t sz = space.getSimpleExtentNpoints();
    assert(sz == 1);
    const DataType &type = H5typeMap<T>::type();
    attribute.read(type, &v);
  }

  //Read a generic array. Array type must be able to provide a pointer to a contiguous memory region.
  //Requires a policy class specifying:
  //ArrayType, ElementType (must be a native type),  static void resize(ArrayType &, const hsize_t sz),  static ElementType* pointer(ArrayType &v)
  template<typename ArrayPolicy, IF_NATIVE(typename ArrayPolicy::ElementType)>
  void read(typename ArrayPolicy::ArrayType &v, const std::string &name){
    typedef typename ArrayPolicy::ElementType T;    
    using namespace H5;
    const DataType &type = H5typeMap<T>::type();
    if(group.back().attrExists(name.c_str())){
      Attribute attribute;
      try{
	 attribute = group.back().openAttribute(name.c_str());
      }catch(H5::Exception& e){
	error_exit(std::cout << "HDF5reader::read(vector)  An attribute of name " << name << " does not exist!\n");
      }  
      DataSpace space = attribute.getSpace();
      hsize_t sz = space.getSimpleExtentNpoints();
      ArrayPolicy::resize(v,sz);
      attribute.read(type, ArrayPolicy::pointer(v));
    }else{
      DataSet dset;
      try{
	dset = group.back().openDataSet(name.c_str());
      }catch(H5::Exception& e){
	error_exit(std::cout << "HDF5reader::read(vector)  A dataset of name " << name << " does not exist!\n");
      }  
      DataSpace space = dset.getSpace();
      hsize_t sz = space.getSimpleExtentNpoints();
      ArrayPolicy::resize(v,sz);
      dset.read(ArrayPolicy::pointer(v), type);
    }
  }

  void enter(const std::string &nm){
    try{
      group.push_back(group.back().openGroup(nm.c_str()));
    }catch(H5::Exception& e){
      error_exit(std::cout << "HDF5reader::enter  A group of name " << nm << " does not exist!\n");
    }
  }
  void leave(){
    group.pop_back();
  }

  inline bool containsGroup(const std::string &name) const{
    try{
      H5::Group grp = group.back().openGroup(name.c_str());
    }catch(H5::Exception& e){
      return false;
    }
    return true;
  }

  bool contains(const std::string &name) const{
    using namespace H5;
    if(group.back().attrExists(name.c_str())){
      return true;
    }else{
      //Check if DataSet
      bool is_dataset = true;
      try{
	DataSet dset = group.back().openDataSet(name.c_str());
      }catch(H5::Exception& e){
	is_dataset = false;
      }
      if(is_dataset) return true;

      //Check if Group
      bool is_group = true;
      try{
	Group grp = group.back().openGroup(name.c_str());
      }catch(H5::Exception& e){
	is_group = false;
      }
      if(is_group) return true;

      return false;
    }
  }
};

CPSFIT_END_NAMESPACE

#endif

#endif
