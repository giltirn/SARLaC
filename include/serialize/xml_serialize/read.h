#ifndef _XML_READ_FUNCS_H___
#define _XML_READ_FUNCS_H___

//Functions to read from an XML file into various object types

#include<serialize/xml_serialize/xml_reader.h>

CPSFIT_START_NAMESPACE

template<typename T, typename std::enable_if<std::is_pod<T>::value, int>::type = 0>
inline void read(XMLreader &reader, T &into, const std::string &tag){
  reader.read(into,tag);
}
inline void read(XMLreader &reader, std::string &into, const std::string &tag){
  reader.read(into,tag);
}
template<typename T, typename std::enable_if<std::is_pod<T>::value, int>::type = 0> //POD version
inline void read(XMLreader &reader, std::vector<T> &v, const std::string &tag){
  reader.read(v,tag);
}
template<typename T, typename std::enable_if<!std::is_pod<T>::value, int>::type = 0> //non POD version
void read(XMLreader &reader, std::vector<T> &v, const std::string &tag){
  reader.enter(tag);
  T tmp;
  for(auto it=reader.group_begin(); it != reader.group_end(); it++){
    if(it->first != "elem") error_exit(std::cout << "read(XMLreader &reader, std::vector<" << printType<T>() << " &v, const std::string &tag) Expected to be in a group full of 'elem' tags. Instead got tag '"<< it->first << "'. Group contains:\n" << reader.printGroupEntries() << std::endl);
    reader.enter(it);
    read(reader,tmp,"");
    v.push_back(tmp);
    reader.leave();
  }
  reader.leave();
}

CPSFIT_END_NAMESPACE

#endif
