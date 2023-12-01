#ifndef _CONSTRAINED_MEMORY_VECTOR_H_
#define _CONSTRAINED_MEMORY_VECTOR_H_

//A vector-like container that if the total memory used by all active instance surpasses some threshold, older vectors are dumped to disk and reloaded later if needed
#include<config.h>
#include<serialize/hdf5_serialize.h>
#include<containers/constrained_memory_vector/memory_manager.h>

SARLAC_START_NAMESPACE

template<typename T>
class constrainedMemoryVector{
  static_assert(is_scalar<T>::value, "Only works for scalar types");
  size_t sz;
  size_t uid;
public:
  constrainedMemoryVector(): sz(0), uid(-1){
  }
  constrainedMemoryVector(const size_t n): constrainedMemoryVector(){
    resize(n);
  }
  constrainedMemoryVector(const size_t n, const T &init): constrainedMemoryVector(n){
    if(n>0){
      T * data = (T *)constrainedMemoryManager::get().get(uid);
      for(int i=0;i<n;i++) data[i] = init;
    }
  }
  constrainedMemoryVector(const std::vector<T> &r): constrainedMemoryVector(r.size()){
    if(r.size() > 0){
      T * data = (T *)constrainedMemoryManager::get().get(uid);
      for(int i=0;i<r.size();i++) data[i] = r[i];
    }
  }
  
  constrainedMemoryVector(const constrainedMemoryVector &r): constrainedMemoryVector(){
    *this = r;
  }
  constrainedMemoryVector(constrainedMemoryVector &&r): sz(r.sz), uid(r.uid){
    r.uid = -1;
  }    

  constrainedMemoryVector & operator=(const constrainedMemoryVector &r){
    this->resize(r.size());
    if(r.size() > 0){
      T * datal = (T *)constrainedMemoryManager::get().get(uid);
      T const* datar = (T const*)constrainedMemoryManager::get().get(r.uid);
      for(int i=0;i<r.size();i++) datal[i] = datar[i];
    }
    return *this;
  } 
  constrainedMemoryVector & operator=(constrainedMemoryVector &&r){
    uid = r.uid;
    sz = r.sz;
    r.uid = -1;   
    return *this;
  } 

  inline const T &operator[](const size_t i) const{
    T const* data = (T const*)constrainedMemoryManager::get().get(uid);
    return data[i];
  }
  inline T &operator[](const size_t i){
    T * data = (T *)constrainedMemoryManager::get().get(uid);
    return data[i];
  }

  inline size_t size() const{ return sz; }
  
  void resize(const size_t _sz){
    if(_sz != sz){
      if(uid != -1) constrainedMemoryManager::get().free(uid);
      uid = constrainedMemoryManager::get().alloc(_sz * sizeof(T));
      //printf("constrainedMemoryVector %p got uid %lu in resize to size %lu\n", this, uid, _sz); fflush(stdout);
      sz = _sz;
    }
  }
  
  ~constrainedMemoryVector(){
    if(uid != -1){
      //printf("constrainedMemoryVector %p with uid %lu and size %lu called free\n", this, uid, sz); fflush(stdout);
      constrainedMemoryManager::get().free(uid);
    }
  }
};

#ifdef HAVE_HDF5
template<typename T>
void write(HDF5writer &writer, const constrainedMemoryVector<T> &value, const std::string &tag){
  writer.enter(tag);
  std::vector<T> tmp(value.size());
  for(int i=0;i<value.size();i++) tmp[i] = value[i];
  write(writer,tmp,"data");
  writer.leave();
}
template<typename T>
void read(HDF5reader &reader, constrainedMemoryVector<T> &value, const std::string &tag){
  reader.enter(tag);
  std::vector<T> tmp;
  read(reader,tmp,"data");
  value.resize(tmp.size());
  for(int i=0;i<tmp.size();i++) value[i] = tmp[i];  
  reader.leave();
}
#endif

SARLAC_END_NAMESPACE

#endif
