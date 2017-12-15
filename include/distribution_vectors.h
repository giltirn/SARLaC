#ifndef _DISTRIBUTION_VECTORS_H_
#define _DISTRIBUTION_VECTORS_H_

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <deque>
#include <unordered_map>
#include <template_wizardry.h>


template<typename T>
using basic_vector = std::vector<T>;

struct _memmappedFileVector_uid{
  inline static size_t & Uid(){ static size_t uid = 0; return uid; }
  inline static size_t getUid(){ return Uid()++; }
};

template<typename T>
class memmappedFileVector{
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
  memmappedFileVector(const basic_vector<T> &r): memmappedFileVector(r.size()){
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



struct mmfile{
  std::string filename;
  boost::iostreams::mapped_file mmp;
  size_t bytes;
  bool created;
  mmfile(): bytes(0), created(false){}
};
class memmappedFileHandleManager{
  std::vector<mmfile*> mmfiles;
  std::deque<int> active;
  const static int max_open = 300;
public:

  //Open the file as a memmapped region
  boost::iostreams::mapped_file & open(const int uid){
    mmfile &mm = *mmfiles[uid];
    assert(mm.created);
    if(!mm.mmp.is_open()){    
#pragma omp critical
      {
	if(active.size() == max_open){
	  mmfiles[active.front()]->mmp.close();
	  active.pop_front();
	}
	mm.mmp.open(mm.filename, boost::iostreams::mapped_file::readwrite, mm.bytes);
	assert(mm.mmp.is_open());
	active.push_back(uid);
      }
    }
    return mm.mmp;
  }

  //Assign a uid and a new entry
  int assign(){
    int uid;
#pragma omp critical
    {
      uid = mmfiles.size();
      mmfiles.resize(mmfiles.size()+1);
      std::ostringstream os; os << "memmappedFileVector_tmp." << uid;
      mmfiles.back() = new mmfile;
      mmfiles.back()->filename = os.str();
      mmfiles.back()->created = false;
    }
    return uid;
  }

  //Create the file on disk
  void setSize(const size_t bytes, const int uid){
    mmfile &mm = *mmfiles[uid];
    if(mm.created){ //resize existing file if necessary
      if(mm.bytes != bytes){
	bool reopen = false;
	if(mm.mmp.is_open()){
	  mm.mmp.close();
	  reopen = true;
	}
	boost::filesystem::resize_file(mm.filename, bytes);
	mm.bytes = bytes;
	if(reopen) mm.mmp.open(mm.filename, boost::iostreams::mapped_file::readwrite, bytes);
	return;
      }
    }else{ //create the file and set its size
      if(!boost::filesystem::exists(mm.filename)){
	std::ofstream of(mm.filename); 
	if(of.fail()){ error_exit(std::cout << "memmappedFileHandleManager failed to open file " << mm.filename << " because " << strerror(errno) << std::endl); }
      }
      if(!boost::filesystem::exists(mm.filename)) error_exit(std::cout << "memmappedFileHandleManager file " << mm.filename << " not open!\n");
      boost::filesystem::resize_file(mm.filename, bytes);
      mm.bytes = bytes;
      mm.created = true;
    }
  }

  inline void closeAndDeleteFile(const int uid){
    mmfile &mm = *mmfiles[uid];
    if(mm.mmp.is_open()) mm.mmp.close();
    boost::filesystem::remove(mm.filename);
    mm.created = false;
    mm.bytes = 0;
  }
    
  ~memmappedFileHandleManager(){
    for(int i=0;i<mmfiles.size();i++){
      closeAndDeleteFile(i);
      delete mmfiles[i];
    }
      
  }

  inline static memmappedFileHandleManager & get(){ static memmappedFileHandleManager mp; return mp; }
};




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
  memmappedManagedFileVector(const basic_vector<T> &r): memmappedManagedFileVector(r.size()){
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


class constrainedMemoryManager{
  struct memblock{
    size_t bytes;
    char* ptr;
    int uid;
    std::string filename;
    size_t last_touch;
    bool ondisk;
    bool delete_file_on_free;

    memblock(): bytes(0), ptr(NULL), filename(""), ondisk(false), delete_file_on_free(false){
#pragma omp critical
      {
	uid = _memmappedFileVector_uid::getUid();
      }
      std::ostringstream os; os << "swapped_mem." << uid;
      filename = os.str();
    }

    memblock(const size_t _bytes): memblock(){ alloc(_bytes); }

    memblock(memblock &&r): bytes(r.bytes), ptr(r.ptr), uid(r.uid), filename(std::move(r.filename)), last_touch(r.last_touch), ondisk(r.ondisk), delete_file_on_free(r.delete_file_on_free){
      r.ptr = NULL;
      r.delete_file_on_free = false;
    }      

    memblock &operator=(memblock &&r){
      assert(ptr == NULL);
      bytes = r.bytes;
      ptr = r.ptr;
      uid = r.uid;
      filename = std::move(r.filename);
      last_touch = r.last_touch;
      ondisk = r.ondisk;
      delete_file_on_free = r.delete_file_on_free;
      r.ptr = NULL;
      r.delete_file_on_free = false;
      return *this;
    }
    
    inline char* alloc(const size_t _bytes){
      //(std::cout << "memblock::alloc bytes " << _bytes << std::endl).flush();
      bytes = _bytes;
      ptr = (char*)malloc(bytes*sizeof(char));
      return ptr;
    }
    
    void write(){
      if(ondisk) return;
      //(std::cout << "memblock::write bytes " << bytes << " to " << filename << std::endl).flush();
      std::ofstream f(filename.c_str(), std::ios::out | std::ios::binary);
      f.write(ptr, bytes);
      assert(!f.bad());
      f.close();
      ::free(ptr); ptr = NULL;
      ondisk = true;
      delete_file_on_free = true;
    }
    void read(){
      //(std::cout << "memblock::read bytes " << bytes << " from " << filename << std::endl).flush();
      assert(ondisk);
      ptr = (char*)malloc(bytes);
      std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
      f.read(ptr, bytes);
      assert(!f.bad());
      f.close();
      ondisk = false;
    }

    inline char* get(){
      static size_t calls = 0;
      #pragma omp critical
      {
	last_touch = calls++;
      }      
      if(ondisk) read();
      return ptr;
    }
    
    void free(){
      if(ptr != NULL){
	///(std::cout << "memblock::free freeing " << bytes << " at ptr " << static_cast<void*>(ptr) << std::endl).flush();
	std::free(ptr);
      }
      if(delete_file_on_free){	
	//(std::cout << "memblock::free deleting file " << filename << std::endl).flush();
	boost::filesystem::remove(filename);
      }
    }
    ~memblock(){
      this->free();
    }      
  };
  typedef std::unordered_map<size_t,memblock> blockMapType;
  blockMapType blocks;
  size_t size;
public:
  
  static inline size_t & maxSize(){
    static size_t max = 10*1024*1024*size_t(1024); //10GB default
    return max;
  }

  void freeUp(const size_t bytes, const size_t skip = -1){
    //(std::cout << "constrainedMemoryManager::freeUp current size " << size << " freeing " << bytes << std::endl).flush();
    size_t bytes_freed = 0;
    typedef std::pair<size_t,blockMapType::iterator> pairType;
    std::list<pairType> its;
    for(blockMapType::iterator it = blocks.begin(); it != blocks.end(); it++) its.push_back(pairType(it->second.last_touch,it));
    its.sort([](const pairType &a, const pairType &b){ return a.first < b.first; });

    for(std::list<pairType>::iterator i = its.begin(); i != its.end(); i++){
      blockMapType::iterator bit = i->second;
      if(bit->second.uid == skip) continue;
      
      bit->second.write();
      size_t fbytes = bit->second.bytes;
      size -= fbytes;
      bytes_freed += fbytes;

      //(std::cout << "constrainedMemoryManager::freeUp freed " << fbytes << " from uid " << bit->second.uid << " increasing total freed to " << bytes_freed << std::endl).flush();
      
      if(bytes_freed >= bytes) return;
    }
    error_exit(std::cout << "Could not free up " << bytes << " bytes!\n");
  }
  
  inline size_t alloc(const size_t bytes){
    //(std::cout << "constrainedMemoryManager::alloc " << bytes << " current size " << size << " max size " << maxSize() << std::endl).flush();
    if(bytes + size > maxSize()) freeUp(bytes + size - maxSize());
    memblock b(bytes);
    size_t uid = b.uid;
    blocks[uid] = std::move(b);
    size += bytes;
    //(std::cout << "constrainedMemoryManager::alloc " << bytes << " assigned uid " << uid << std::endl).flush();
    return uid;
  }

  inline char* get(const size_t uid){
    blockMapType::iterator it = blocks.find(uid);
    if(it == blocks.end()) error_exit(std::cout << "constrainedMemoryManager::get Could not find uid " << uid << " (present largest uid=" << _memmappedFileVector_uid::Uid() << ")\n");

    //(std::cout << "constrainedMemoryManager::get uid " << uid << std::endl).flush();
    if(it->second.ondisk && size + it->second.bytes > maxSize()){
      freeUp(size + it->second.bytes - maxSize(), uid);
      size += it->second.bytes;
    }
    return it->second.get();
  }
  
  void free(const size_t uid){
    blockMapType::iterator it = blocks.find(uid);
    if(it == blocks.end()) error_exit(std::cout << "constrainedMemoryManager::free Could not find uid " << uid << " (present largest uid=" << _memmappedFileVector_uid::Uid() << ")\n");
    
    size -= it->second.bytes;
    blocks.erase(it);	
  }

  inline static constrainedMemoryManager & get(){ static constrainedMemoryManager mp; return mp; }
};



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
  constrainedMemoryVector(const basic_vector<T> &r): constrainedMemoryVector(r.size()){
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


#endif
