#ifndef _MEMMAPPED_FILE_VECTOR_H_
#define _MEMMAPPED_FILE_VECTOR_H_

//A container for scalar types that maps directly to a file on disk

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <config.h>
#include <utils/macros.h>
#include <serialize/hdf5_serialize.h>

CPSFIT_START_NAMESPACE

template<typename T>
class memmappedFileVector{
  struct _memmappedFileVector_uid{
    inline static size_t & Uid(){ static size_t uid = 0; return uid; }
    inline static size_t getUid(){ return Uid()++; }
  };

  static_assert(is_scalar<T>::value, "Only works for scalar types");
  std::string filename;
  boost::iostreams::mapped_file file;
  size_t bytes;
  size_t sz;
public:
  memmappedFileVector(): sz(0){
    int uid;
#pragma omp critical
    {
      uid = _memmappedFileVector_uid::getUid();
    }
    std::ostringstream os; os << "memmappedFileVector_tmp." << uid;
    filename = os.str();
  }
  memmappedFileVector(const size_t n): memmappedFileVector(){
    resize(n);
  }
  memmappedFileVector(const size_t n, const T &init): memmappedFileVector(n){
    T * data = (T *)file.data();
    for(int i=0;i<n;i++) data[i] = init;
  }
  memmappedFileVector(const std::vector<T> &r): memmappedFileVector(r.size()){
    for(int i=0;i<r.size();i++) this->operator[](i) = r[i];
  }
  
  memmappedFileVector(const memmappedFileVector &r): memmappedFileVector(r.size()){
    for(int i=0;i<r.size();i++) this->operator[](i) = r[i];
  }
  memmappedFileVector(memmappedFileVector &&r): filename(std::move(r.filename)), bytes(r.bytes), sz(r.sz){
    r.file.close();
    file.open(filename, boost::iostreams::mapped_file::readwrite, bytes);
    assert(file.is_open());
  }    

  memmappedFileVector & operator=(const memmappedFileVector &r){ 
    this->resize(r.size());
    for(int i=0;i<r.size();i++) this->operator[](i) = r[i];
    return *this;
  } 
  memmappedFileVector & operator=(memmappedFileVector &&r){ 
    r.file.close();
    filename = std::move(r.filename);
    bytes=r.bytes;
    sz = r.sz;
    file.open(filename, boost::iostreams::mapped_file::readwrite, bytes);
    assert(file.is_open());    
    return *this;
  } 

  inline const T &operator[](const size_t i) const{ return *((T *)file.const_data()+i); }
  inline T &operator[](const size_t i){ return *((T*)file.data()+i); }

  inline size_t size() const{ return sz; }
  
  void resize(const size_t _sz){
    if(_sz == sz) return;
    
    if(file.is_open()) file.close();    
    sz = _sz;
    bytes = _sz * sizeof(T);
    
    if(!boost::filesystem::exists(filename)){
      std::ofstream of(filename); 
      if(of.fail()){
	error_exit(std::cout << "memmappedFileVector failed to open file " << filename << " because " << strerror(errno) << std::endl);
      }
    }
    if(!boost::filesystem::exists(filename)) error_exit(std::cout << "memmappedFileVector file " << filename << " not open!\n");
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

//HDF5 serialization routines
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


CPSFIT_END_NAMESPACE

#endif
