#ifndef _MEMMAPPED_MANAGED_FILE_VECTOR_HANDLE_MANAGER_H__
#define _MEMMAPPED_MANAGED_FILE_VECTOR_HANDLE_MANAGER_H__

#include<sstream>
#include<iostream>

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <config.h>
#include <utils/macros.h>
#include <utils/utils.h>

SARLAC_START_NAMESPACE

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

SARLAC_END_NAMESPACE

#endif
