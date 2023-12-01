#ifndef _SARLAC_QDP_BINARY_READ_H__
#define _SARLAC_QDP_BINARY_READ_H__

//Allow reading of files written in QDP's binary format, eg from UKvalence
#include<fstream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils/endian.h>
#include<utils/utils/error.h>

SARLAC_START_NAMESPACE

class QDPbinaryReader{
  bool system_bigendian;
  std::ifstream strm;
  std::string file;
public:
  QDPbinaryReader(const std::string &file): file(file), system_bigendian(isBigEndian()), strm(file.c_str(), std::ios::binary){    
    if(!strm.good()) error_exit(std::cout << "QDPbinaryReader could not open file " << file << std::endl);
    strm.exceptions ( std::ios::failbit | std::ios::badbit );
  }
  
  template<typename T>
  typename std::enable_if<std::is_fundamental<T>::value, void>::type read(T &into){
    strm.read((char*)&into, sizeof(T));
    if(!system_bigendian)
      into = reverseEndianness(into);
  }

  void read(std::string &into){
    std::getline(strm,into);
  }

  const std::string &fileName() const{ return file; }
};

template<typename T>
inline typename std::enable_if<std::is_fundamental<T>::value, void>::type read(QDPbinaryReader &rd, T &into){
  rd.read(into);
}

inline void read(QDPbinaryReader &rd, std::string &into){
  rd.read(into);
}

template<typename T>
inline void read(QDPbinaryReader &rd, std::vector<T> &into){
  int size;
  rd.read(size);

  into.resize(size);
  for(int i=0;i<size;i++) read(rd,into[i]);
}

SARLAC_END_NAMESPACE


#endif
