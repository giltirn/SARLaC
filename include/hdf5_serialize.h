#ifndef _HDF5_SERIALIZE_H_
#define _HDF5_SERIALIZE_H_

#include<config.h>

#ifdef HAVE_HDF5
#include <map>
#include <H5Cpp.h>
#include <boost/preprocessor.hpp>

#include<utils.h>



template<typename T>
struct H5typeMap{
  enum {is_native = 0};
};

#define H5_TYPE_MAP(TYPE, TYPENAME)					\
  template<>								\
  struct H5typeMap<TYPE>{						\
    static inline const H5::DataType & type(void){ return H5::PredType::TYPENAME; } \
    enum {is_native = 1};						\
  }

H5_TYPE_MAP(bool, NATIVE_B8);
H5_TYPE_MAP(char, NATIVE_CHAR);
H5_TYPE_MAP(signed char, NATIVE_SCHAR);
H5_TYPE_MAP(unsigned char, NATIVE_UCHAR);
H5_TYPE_MAP(short, NATIVE_SHORT);
H5_TYPE_MAP(unsigned short, NATIVE_USHORT);
H5_TYPE_MAP(int, NATIVE_INT);
H5_TYPE_MAP(unsigned int, NATIVE_UINT);
H5_TYPE_MAP(long, NATIVE_LONG);
H5_TYPE_MAP(unsigned long, NATIVE_ULONG);
H5_TYPE_MAP(long long, NATIVE_LLONG);
H5_TYPE_MAP(unsigned long long, NATIVE_ULLONG);
H5_TYPE_MAP(float, NATIVE_FLOAT);
H5_TYPE_MAP(double, NATIVE_DOUBLE);
H5_TYPE_MAP(long double, NATIVE_LDOUBLE);


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
  
class HDF5reader{
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

public:
  HDF5reader(const std::string &filename): file(filename.c_str(), H5F_ACC_RDONLY){
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
  void read(std::string &v, const std::string &name){
    read<HDF5readerStringPolicy>(v,name);
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


//Strings
inline void write(HDF5writer &writer, const std::string &value, const std::string &tag){
  writer.write(value,tag);
}
inline void read(HDF5reader &reader, std::string &value, const std::string &tag){
  reader.read(value,tag);
}



template<typename T, int>
struct _isDistributionOfHDF5nativetype{ enum {value = 0}; };

template<typename T>
struct _isDistributionOfHDF5nativetype<T,1>{ enum {value = H5typeMap<typename T::DataType>::is_native}; }; 

template<typename T>
struct isDistributionOfHDF5nativetype{ enum {value = _isDistributionOfHDF5nativetype<T,hasSampleMethod<T>::value>::value }; };    

  
//Non-native vectors (but not distributions, they are specialized)
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

//Pair
template<typename T, typename U>
inline static void write(HDF5writer &writer, const std::pair<T,U> &value, const std::string &tag){
  writer.enter(tag); //enter a group
  write(writer, value.first, "first");
  write(writer, value.second, "second");
  writer.leave();
}
template<typename T, typename U>
inline static void read(HDF5reader &reader, std::pair<T,U> &value, const std::string &tag){
  reader.enter(tag); //enter a group
  read(reader, value.first, "first");
  read(reader, value.second, "second");
  reader.leave();
}

//Map
template<typename KeyT, typename DataT>
static void write(HDF5writer &writer, const std::map<KeyT,DataT> &value, const std::string &tag){
  writer.enter(tag);
  int sz = value.size();
  write(writer,sz,"sz");
  writer.enter("entries");
  int i=0;
  for(typename std::map<KeyT,DataT>::const_iterator it = value.begin(); it != value.end(); it++, i++){
    const std::pair<KeyT,DataT> &e = *it;
    std::ostringstream os; os << "entry_" << i;
    write(writer,e,os.str());
  }
  writer.leave();
  writer.leave();
}
template<typename KeyT, typename DataT>
static void read(HDF5reader &reader, std::map<KeyT,DataT> &value, const std::string &tag){
  value.clear();
  reader.enter(tag);
  int sz;
  read(reader,sz,"sz");
  reader.enter("entries");
  for(int i=0;i<sz;i++){
    std::pair<KeyT,DataT> e;
    std::ostringstream os; os << "entry_" << i;
    read(reader,e,os.str());
    value.emplace(std::move(e));
  }
  assert(value.size() == sz);
  reader.leave();
  reader.leave();
}


#define _GENERATE_HDF5_SERIALIZE_METHOD_WRITEIT(r,data,elem) ::write(writer, this-> elem, BOOST_PP_STRINGIZE(elem)); 
#define _GENERATE_HDF5_SERIALIZE_METHOD_READIT(r,data,elem) ::read(reader, this-> elem, BOOST_PP_STRINGIZE(elem)); 


//Generate HDF5 serialize methods for a class. MEMBERS should be a series of member names in parentheses, eg (member1)(member2)(member3)....
#define GENERATE_HDF5_SERIALIZE_METHOD(MEMBERS)\
    void write(HDF5writer &writer, const std::string &tag) const{ \
      writer.enter(tag);\
      BOOST_PP_SEQ_FOR_EACH(_GENERATE_HDF5_SERIALIZE_METHOD_WRITEIT, , MEMBERS); \
      writer.leave(); \
    }\
    void read(HDF5reader &reader, const std::string &tag){ \
      reader.enter(tag); \
      BOOST_PP_SEQ_FOR_EACH(_GENERATE_HDF5_SERIALIZE_METHOD_READIT, , MEMBERS); \
      reader.leave(); \
    }

#define GENERATE_HDF5_SERIALIZE_FUNC(CLASSNAME)\
  inline void write(HDF5writer &writer, const CLASSNAME &d, const std::string &tag){ d.write(writer,tag); }\
  inline void read(HDF5reader &reader, CLASSNAME &d, const std::string &tag){ d.read(reader,tag); }



//Example for a struct/class
/*
struct S{
  int a;
  int b;
  int c;

  S(const int _a, const int _b, const int _c): a(_a), b(_b), c(_c){}
};

void write(HDF5writer &writer, const S &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.a,"a");
  write(writer,value.b,"b");
  write(writer,value.c,"c");
  writer.leave();
}
void read(HDF5reader &reader, S &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.a,"a");
  read(reader,value.b,"b");
  read(reader,value.c,"c");
  reader.leave();
}
*/


#else

//Empty macros so we don't need to keep adding #ifdef HAVE_HDF5 around everything
#define GENERATE_HDF5_SERIALIZE_METHOD(MEMBERS)
#define GENERATE_HDF5_SERIALIZE_FUNC(CLASSNAME)


#endif

#endif
