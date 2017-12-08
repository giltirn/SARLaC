#ifndef _DISTRIBUTION_VECTORS_H_
#define _DISTRIBUTION_VECTORS_H_

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <template_wizardry.h>


template<typename T>
using basic_vector = std::vector<T>;

struct _memmappedFileVector_uid{
  inline static int getUid(){ static int uid = 0; return uid++; }
};

template<typename T>
class memmappedFileVector{
  static_assert(is_scalar<T>::value, "Only works for scalar types");
  std::string filename;
  boost::iostreams::mapped_file file;
  size_t bytes;
  size_t sz;
public:
  memmappedFileVector(){
    std::ostringstream os; os << "memmappedFileVector_tmp." << _memmappedFileVector_uid::getUid();
    filename = os.str();
  }
  memmappedFileVector(const size_t n): memmappedFileVector(){
    resize(n);
  }
  memmappedFileVector(const size_t n, const T &init): memmappedFileVector(n){
    T * data = (T *)file.data();
    for(int i=0;i<n;i++) data[i] = init;
  }
  memmappedFileVector(const basic_vector<T> &r): memmappedFileVector(r.size()){
    for(int i=0;i<r.size();i++) this->operator[](i) = r[i];
  }
    
  
  inline const T &operator[](const size_t i) const{ return *((T *)file.const_data()+i); }
  inline T &operator[](const size_t i){ return *((T*)file.data()+i); }

  inline size_t size() const{ return sz; }
  
  void resize(const size_t _sz){
    if(_sz == sz) return;
    
    if(file.is_open()) file.close();    
    sz = _sz;
    bytes = _sz * sizeof(T);
    
    { std::ofstream of(filename); }
    boost::filesystem::resize_file(filename, bytes);
    
    file.open(filename, boost::iostreams::mapped_file::readwrite, bytes);
    assert(file.is_open());
  }
  
  ~memmappedFileVector(){
    if(file.is_open()){
      file.close();
      boost::filesystem::remove(filename);
    } 
  }
};

#ifdef HAVE_HDF5
template<typename T>
void write(HDF5writer &writer, const memmappedFileVector<T> &value, const std::string &tag){
  writer.enter(tag);
  std::vector<T> tmp(value.size());
  for(int i=0;i<value.size();i++) tmp[i] = value[i];
  write(writer,tmp,"data");
  writer.leave();
}
template<typename T>
void read(HDF5reader &reader, memmappedFileVector<T> &value, const std::string &tag){
  reader.enter(tag);
  std::vector<T> tmp;
  read(reader,tmp,"data");
  value.resize(tmp.size());
  for(int i=0;i<tmp.size();i++) value[i] = tmp[i];  
  reader.leave();
}
#endif

#endif
