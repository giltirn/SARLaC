
#ifndef _CPSFIT_TAGGED_VALUE_CONTAINER_H_
#define _CPSFIT_TAGGED_VALUE_CONTAINER_H_

//A vector-like container with all the usual fit boilerplate that has a shared tag->index mapping allowing for, eg. dynamic parameter containers with named parameters

#include<unordered_map>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<ET/generic_ET.h>

CPSFIT_START_NAMESPACE

template<typename T, typename TagType = std::string>
class taggedValueContainer{
  std::unordered_map<TagType,size_t> const* tag_map; //unique and shared over instances of this class
  std::vector<T> t;
public:
  ENABLE_GENERIC_ET(taggedValueContainer, taggedValueContainer<T>, taggedValueContainer<T>);
  taggedValueContainer(): tag_map(NULL), t(0){}
  taggedValueContainer(std::unordered_map<TagType,size_t> const* tag_map): t(tag_map->size()), tag_map(tag_map){}

  //Integer accessors
  inline T & operator()(const size_t i){ return t[i]; }
  inline const T &operator()(const size_t i) const{ return t[i]; }

  //Tagged accessors
  inline T & operator()(const TagType &i){ assert(tag_map != NULL); auto it = tag_map->find(i); assert(it != tag_map->end()); return t[it->second]; }
  inline const T &operator()(const TagType &i) const{ assert(tag_map != NULL); auto it = tag_map->find(i); assert(it != tag_map->end()); return t[it->second]; }
  
  //This version also returns the index
  inline T & operator()(size_t &idx, const TagType &i){ assert(tag_map != NULL); auto it = tag_map->find(i); assert(it != tag_map->end()); idx = it->second; return t[idx]; }
  inline const T &operator()(size_t &idx, const TagType &i) const{ assert(tag_map != NULL); auto it = tag_map->find(i); assert(it != tag_map->end()); idx = it->second; return t[idx]; }
  
  //Usual boilerplate stuff
  inline size_t size() const{ return t.size(); }
  inline std::string print() const{ 
    if(tag_map == NULL) return "";
    std::vector<typename std::unordered_map<TagType,size_t>::const_iterator> inv_map(tag_map->size());
    for(auto it = tag_map->begin(); it != tag_map->end(); it++) inv_map[it->second] = it;
    std::ostringstream os;  os << "{";
    for(int i=0;i<tag_map->size();i++)
      os << (i > 0 ? ", " : "") << inv_map[i]->first << "=" << t[i];
    os << "}";    
    return os.str(); 
  }
  inline void resize(const size_t i){ if(i!=t.size()) error_exit(std::cout << printType<taggedValueContainer<T> >() << " resize called with value " << i << " != " << t.size() << "\n"); }

  //This version is used by the ETE but can also be used manually to setup an instance generated with the default constructor
  inline void resize(std::unordered_map<TagType,size_t> const* mp){ tag_map = mp; t.resize(mp->size()); }
  
  inline void zero(){ for(auto it=t.begin(); it != t.end(); ++it) *it=0.; }

  //Access the tag map
  inline std::unordered_map<TagType,size_t> const* getTagMap() const{ return tag_map; }
  
  //Parameter index
  size_t index(const TagType &tag) const{ auto it = tag_map->find(tag); assert(it != tag_map.end()); return it->second; }
};

template<typename T, typename TagType>
struct getElem<taggedValueContainer<T,TagType> >{
  static inline T& elem(taggedValueContainer<T,TagType> &v, const int i){ return v(i); }
  static inline const T& elem(const taggedValueContainer<T,TagType> &v, const int i){ return v(i); }    
  inline static std::unordered_map<TagType,size_t> const* common_properties(const taggedValueContainer<T,TagType> &v){ return v.getTagMap(); } //for setting up ET output
};
template<typename T, typename TagType>
std::ostream & operator<<(std::ostream &os, const taggedValueContainer<T,TagType> &r){
  os << r.print(); return os;
}

CPSFIT_END_NAMESPACE

#endif
