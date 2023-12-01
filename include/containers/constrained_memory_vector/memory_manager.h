#ifndef _CONSTRAINED_MEMORY_MANAGER_H_
#define _CONSTRAINED_MEMORY_MANAGER_H_

#include<iostream>
#include<sstream>
#include<list>
#include<fstream>
#include<unordered_map>
#include<cassert>

#include <boost/filesystem.hpp>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>

SARLAC_START_NAMESPACE


class constrainedMemoryManager{
  struct memblock{
    size_t bytes;
    char* ptr;
    int uid;
    std::string filename;
    size_t last_touch;
    bool ondisk;
    bool delete_file_on_free;

    inline static size_t & Uid(){ static size_t uid = 0; return uid; }
    inline static size_t getUid(){ return Uid()++; }

    memblock(): bytes(0), ptr(NULL), filename(""), ondisk(false), delete_file_on_free(false){
#pragma omp critical
      {
	uid = getUid();
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
    if(it == blocks.end()) error_exit(std::cout << "constrainedMemoryManager::get Could not find uid " << uid << " (present largest uid=" << memblock::Uid() << ")\n");

    //(std::cout << "constrainedMemoryManager::get uid " << uid << std::endl).flush();
    if(it->second.ondisk && size + it->second.bytes > maxSize()){
      freeUp(size + it->second.bytes - maxSize(), uid);
      size += it->second.bytes;
    }
    return it->second.get();
  }
  
  void free(const size_t uid){
    blockMapType::iterator it = blocks.find(uid);
    if(it == blocks.end()) error_exit(std::cout << "constrainedMemoryManager::free Could not find uid " << uid << " (present largest uid=" << memblock::Uid() << ")\n");
    
    size -= it->second.bytes;
    blocks.erase(it);	
  }

  inline static constrainedMemoryManager & get(){ static constrainedMemoryManager mp; return mp; }
};

SARLAC_END_NAMESPACE

#endif
