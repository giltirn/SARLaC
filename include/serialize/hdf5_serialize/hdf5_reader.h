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
  template<typename T>
  struct HDF5readerVectorPolicy{
    typedef std::vector<T> ArrayType;
    typedef T ElementType;
    inline static void resize(ArrayType &v, const int sz){ v.resize(sz); }
    inline static ElementType* pointer(ArrayType &v){ return v.data(); }
  };
  template<typename T, std::size_t Size>
  struct HDF5readerArrayPolicy{
    typedef std::array<T,Size> ArrayType;
    typedef T ElementType;
    inline static void resize(ArrayType &v, const int sz){ assert(sz == Size); }
    inline static ElementType* pointer(ArrayType &v){ return v.data(); }
  };
  struct HDF5readerStringPolicy{
    typedef std::string ArrayType;
    typedef char ElementType;
    inline static void resize(ArrayType &v, const int sz){ v.resize(sz); }
    inline static ElementType* pointer(ArrayType &v){ return &v[0]; }
  };
  template<typename T>
  struct HDF5readerComplexPolicy{
    typedef std::complex<T> ArrayType;
    typedef T ElementType;
    inline static void resize(ArrayType &v, const int sz){ assert(sz == 2); }
    inline static ElementType* pointer(ArrayType &v){ return reinterpret_cast<T*>(&v); }
  };
  H5::H5File file;
  std::vector<H5::Group> group;

  template<typename ArrayPolicy>
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
  template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
  void read(T &v, const std::string &name){
    using namespace H5;
    Attribute attribute = group.back().openAttribute(name.c_str());
    DataSpace space = attribute.getSpace();
    hsize_t sz = space.getSimpleExtentNpoints();
    assert(sz == 1);
    const DataType &type = H5typeMap<T>::type();
    attribute.read(type, &v);
  }
  template<typename T, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
  inline void read(std::vector<T> &v, const std::string &name){
    read<HDF5readerVectorPolicy<T> >(v,name);
  }
  template<typename T, std::size_t Size, typename std::enable_if<H5typeMap<T>::is_native, int>::type = 0>
  void read(std::array<T,Size> &v, const std::string &name){
    read<HDF5readerArrayPolicy<T,Size> >(v,name);
  }  
  inline void read(std::string &v, const std::string &name){
    read<HDF5readerStringPolicy>(v,name);
  }
  template<typename T>
  inline void read(std::complex<T> &v, const std::string &name){
    read<HDF5readerComplexPolicy<T> >(v,name);
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
