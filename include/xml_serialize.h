#ifndef _XML_SERIALIZE_H_
#define _XML_SERIALIZE_H_

//An XML IO class and routines that has the same functionality as the HDF5 IO
#include <iostream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <utils.h>

//Space-separated POD array
template<typename PODtype, typename std::enable_if<std::is_pod<PODtype>::value, int>::type = 0>
struct pod_vector_translator
{
  typedef std::string internal_type;
  typedef std::vector<PODtype> external_type;
  boost::optional<external_type> get_value(const internal_type& str){
    boost::optional<external_type> out(true,external_type());
    std::stringstream ss(str);
    PODtype tmp;
    while(ss){
      if(ss.eof()) break;
      ss >> tmp;
      out.get().push_back(tmp);
    }
    return out;
  }
};

 
class XMLreader{
  boost::property_tree::ptree pt;
  std::vector<std::string> path;

  std::vector<boost::property_tree::ptree const*> group;

  template<typename T>
  static inline T getEntry(const std::string &path, const boost::property_tree::ptree &from){
    using boost::property_tree::ptree;
    typedef ptree::path_type path_type;
    boost::optional<T> e = from.get_optional<T>(path_type(path));
    if(!e) error_exit(std::cout << "getEntry could not parse type " << printType<T>() << " from path " << path << std::endl);
    return e.get();
  }
  template<typename T,typename Translator>
  static inline T getEntryCustom(const std::string &path, Translator tr, const boost::property_tree::ptree &from){
    using boost::property_tree::ptree;
    typedef ptree::path_type path_type;
    boost::optional<T> e = from.get_optional<T>(path_type(path),tr);
    if(!e) error_exit(std::cout << "getEntryCustom could not parse type " << printType<T>() << " from path \"" << path << "\". Environment contains tags:\n" << printGroupEntries(from) << std::endl);
    return e.get();
  }
  static inline auto getChild(const std::string &path, const boost::property_tree::ptree &from)->decltype( from.get_child_optional(boost::property_tree::ptree::path_type(path)).get() ){
    using boost::property_tree::ptree;
    typedef ptree::path_type path_type;
    auto e = from.get_child_optional(path_type(path));
    if(!e) error_exit(std::cout << "getChild could not parse child node from path " << path << std::endl);
    return e.get();
  }
public:
  XMLreader(const std::string &filename){
    std::ifstream is(filename);
    if(!is.good()) error_exit(std::cout << "XMLreader::XMLreader issue with opening file " << filename << std::endl);
    read_xml(is, pt);
    group.push_back(&pt);
  }

  //Read a single value of POD types
  template<typename T, typename std::enable_if<std::is_pod<T>::value, int>::type = 0>
  void read(T &v, const std::string &name){
    v = getEntry<T>(name, *group.back());
  }
  //Read a string
  void read(std::string &v, const std::string &name){
    v = getEntry<std::string>(name, *group.back());
  }
  
  template<typename T, typename std::enable_if<std::is_pod<T>::value, int>::type = 0>
  void read(std::vector<T> &v, const std::string &name){
    //Vectors of plain-old data. Expect space separated
    pod_vector_translator<T> tr;
    v = getEntryCustom<std::vector<T> >(name, tr, *group.back());
  }
  
  void enter(const std::string &nm){
    if(nm == ""){ group.push_back(group.back()); return; }//just copy existing pointer for empty push request. 
    
    using boost::property_tree::ptree;
    typedef ptree::path_type path_type;
    auto e = group.back()->get_child_optional(path_type(nm));
    if(!e) error_exit(std::cout << "getChild could not parse child node from path \"" << nm << "\". Environment contains tags:\n" << printGroupEntries() << std::endl);
    group.push_back(&e.get());
  }
  // //Version of the above that enters the i'th child with name 'nm'
  // void enter(const std::string &nm, const int i){
  //   int ii=0;
  //   bool found = false;
  //   for(auto it = group.back()->begin(); it != group.back()->end(); it++){
  //     if(it->first == nm){
  // 	if(ii != i) ++ii;
  // 	else{
  // 	  found = true;
  // 	  group.push_back(&it->second);
  // 	  return;
  // 	}
  //     }
  //   }
  //   if(!found) error_exit(std::cout << "XMLreader(" << nm << ", " << i << ") could not find an entry with that index\n");
  // }
  void enter(typename boost::property_tree::ptree::const_iterator it){
    group.push_back(&it->second);
  }
	     
  void leave(){
    group.pop_back();
  }

  static std::string printGroupEntries(const boost::property_tree::ptree &pt){
    std::ostringstream os;
    for(auto it = pt.begin(); it != pt.end(); it++)
      os << it->first << std::endl;
    return os.str();
  }
  std::string printGroupEntries(){ return printGroupEntries(*group.back()); }
  
  typename boost::property_tree::ptree::const_iterator group_begin() const{ return group.back()->begin(); }
  typename boost::property_tree::ptree::const_iterator group_end() const{ return group.back()->end(); }
    
};

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

#endif
