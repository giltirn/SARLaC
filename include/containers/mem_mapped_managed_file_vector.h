#ifndef _MEMMAPPED_MANAGED_FILE_VECTOR_H_
#define _MEMMAPPED_MANAGED_FILE_VECTOR_H_

//A variant of memmappedFileVector that only allows a maximum number to be open at any given time, after which earlier created files are closed (reopened if necessary later)

#include<serialize/hdf5_serialize.h>
#include<containers/mem_mapped_managed_file_vector/handle_manager.h>

CPSFIT_START_NAMESPACE

template<typename T>
class memmappedManagedFileVector{
  static_assert(is_scalar<T>::value, "Only works for scalar types");
  size_t sz;
  int uid;
public:
  memmappedManagedFileVector(): sz(0){
    uid = memmappedFileHandleManager::get().assign(); //assign a uid and filename for this object
  }
  memmappedManagedFileVector(const size_t n): memmappedManagedFileVector(){
    resize(n);
  }
  memmappedManagedFileVector(const size_t n, const T &init): memmappedManagedFileVector(n){
    if(n>0){
      T * data = (T *)memmappedFileHandleManager::get().open(uid).data();
      for(int i=0;i<n;i++) data[i] = init;
    }
  }
  memmappedManagedFileVector(const std::vector<T> &r): memmappedManagedFileVector(r.size()){
    if(r.size() > 0){
      T * data = (T *)memmappedFileHandleManager::get().open(uid).data();
      for(int i=0;i<r.size();i++) data[i] = r[i];
    }
  }
  
  memmappedManagedFileVector(const memmappedManagedFileVector &r): memmappedManagedFileVector(r.size()){
    if(r.size() > 0){
      T * data = (T *)memmappedFileHandleManager::get().open(uid).data();
      for(int i=0;i<r.size();i++) data[i] = r[i];
    }
  }
  memmappedManagedFileVector(memmappedManagedFileVector &&r): sz(r.sz), uid(r.uid){
    r.uid = -1;
  }    

  memmappedManagedFileVector & operator=(const memmappedManagedFileVector &r){
    this->resize(r.size());
    if(r.size() > 0){
      T * data = (T *)memmappedFileHandleManager::get().open(uid).data();
      for(int i=0;i<r.size();i++) data[i] = r[i];
    }
    return *this;
  } 
  memmappedManagedFileVector & operator=(memmappedManagedFileVector &&r){
    uid = r.uid;
    sz = r.sz;
    r.uid = -1;   
    return *this;
  } 

  inline const T &operator[](const size_t i) const{
    T const* data = (T *)memmappedFileHandleManager::get().open(uid).const_data();
    return data[i];
  }
  inline T &operator[](const size_t i){
    T * data = (T *)memmappedFileHandleManager::get().open(uid).data();
    return data[i];
  }

  inline size_t size() const{ return sz; }
  
  void resize(const size_t _sz){
    if(_sz != sz){
      memmappedFileHandleManager::get().setSize(_sz * sizeof(T), uid);
      sz = _sz;
    }
  }
  
  ~memmappedManagedFileVector(){
    if(uid != -1) memmappedFileHandleManager::get().closeAndDeleteFile(uid);
  }
};

#ifdef HAVE_HDF5
template<typename T>
void write(HDF5writer &writer, const memmappedManagedFileVector<T> &value, const std::string &tag){
  writer.enter(tag);
  std::vector<T> tmp(value.size());
  for(int i=0;i<value.size();i++) tmp[i] = value[i];
  write(writer,tmp,"data");
  writer.leave();
}
template<typename T>
void read(HDF5reader &reader, memmappedManagedFileVector<T> &value, const std::string &tag){
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
