#ifndef _CPSFIT_UTILS_FILESYSTEM_H_
#define _CPSFIT_UTILS_FILESYSTEM_H_

#include<fstream>

#include<config.h>
#include<utils/macros.h>


CPSFIT_START_NAMESPACE

inline bool fileExists(const std::string &filename){
  std::ifstream infile(filename);
  return infile.good();
}

CPSFIT_END_NAMESPACE
#endif
