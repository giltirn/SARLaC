#ifndef _SARLAC_UTILS_FILESYSTEM_H_
#define _SARLAC_UTILS_FILESYSTEM_H_

#include<vector>
#include<string>
#include<fstream>
#include<regex>
#include<boost/filesystem.hpp>


#include<config.h>
#include<utils/macros.h>
#include<utils/utils/error.h>

SARLAC_START_NAMESPACE

inline bool fileExists(const std::string &filename){
  std::ifstream infile(filename);
  return infile.good();
}

inline std::vector<std::string> listFiles(const std::string &dir, const std::string &regex_fmt){
  using namespace boost::filesystem;
  path p(dir);
  if(!exists(p)) error_exit(std::cout << "Directory \"" << dir << "\" does not exist!\n");
  if(!is_directory(p)) error_exit(std::cout << "\"" << dir << "\" is not a directory!\n");
  
  std::vector<std::string> out;

  std::regex r(regex_fmt);

  for(directory_entry& x : directory_iterator(p)){
    if(is_regular_file(x.path())){
      std::string file = x.path().stem().string();
      //std::cout << file << std::endl;
      if(std::regex_search(file, r)) out.push_back(file);
    }
  }
  return out;
}

//returns true if a new directory was created
inline bool createDirectory(const std::string &dir){
  return boost::filesystem::create_directories(dir);
}

SARLAC_END_NAMESPACE
#endif
