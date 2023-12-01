
#ifndef _SARLAC_TAGGED_VALUE_CONTAINER_H_
#define _SARLAC_TAGGED_VALUE_CONTAINER_H_

//A vector-like container with all the usual fit boilerplate that has a shared tag->index mapping allowing for, eg. dynamic parameter containers with named parameters
#include<list>
#include<unordered_map>
#include<map>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<ET/generic_ET.h>

SARLAC_START_NAMESPACE

template<typename T, typename TagType = std::string>
class taggedValueContainer{
public:
  typedef std::unordered_map<TagType,size_t> TagIdxMap;
private:
  TagIdxMap const* tag_map; //unique and shared over instances of this class
  std::vector<T> t;

  inline size_t checkGetIdx(const TagType &tag) const{ 
    if(tag_map == NULL) error_exit(std::cout << "taggedValueContainer::checkGetIdx tag map is NULL" << std::endl);
    auto it = tag_map->find(tag); 
    if(it == tag_map->end()) error_exit(std::cout << "taggedValueContainer::checkGetIdx could not find tag " << tag << std::endl);
    return it->second;
  }

public:
  ENABLE_GENERIC_ET(taggedValueContainer, taggedValueContainer<T>, taggedValueContainer<T>);
  taggedValueContainer(): tag_map(NULL), t(0){}
  taggedValueContainer(TagIdxMap const* tag_map): t(tag_map->size()), tag_map(tag_map){}

  //Integer accessors
  inline T & operator()(const size_t i){ return t[i]; }
  inline const T &operator()(const size_t i) const{ return t[i]; }

  //Tagged accessors
  inline T & operator()(const TagType &i){ return t[checkGetIdx(i)]; }
  inline const T &operator()(const TagType &i) const{ return t[checkGetIdx(i)]; }
  
  //This version also returns the index
  inline T & operator()(size_t &idx, const TagType &i){ idx = checkGetIdx(i); return t[idx]; }
  inline const T &operator()(size_t &idx, const TagType &i) const{ idx = checkGetIdx(i); return t[idx]; }
  
  //Usual boilerplate stuff
  inline size_t size() const{ return t.size(); }
  inline std::string print() const{ 
    if(tag_map == NULL) return "";
    std::vector<typename TagIdxMap::const_iterator> inv_map(tag_map->size());
    for(auto it = tag_map->begin(); it != tag_map->end(); it++) inv_map[it->second] = it;
    std::ostringstream os;  os << "{";
    for(int i=0;i<tag_map->size();i++)
      os << (i > 0 ? ", " : "") << inv_map[i]->first << "=" << t[i];
    os << "}";    
    return os.str(); 
  }
  inline void resize(const size_t i){ if(i!=t.size()) error_exit(std::cout << printType<taggedValueContainer<T> >() << " resize called with value " << i << " != " << t.size() << "\n"); }

  //This version is used by the ETE but can also be used manually to setup an instance generated with the default constructor
  inline void resize(TagIdxMap const* mp){ tag_map = mp; t.resize(mp->size()); }
  
  inline void zero(){ for(auto it=t.begin(); it != t.end(); ++it) *it=0.; }

  //Access the tag map
  inline TagIdxMap const* getTagMap() const{ return tag_map; }
  
  //Parameter index
  size_t index(const TagType &tag) const{ return checkGetIdx(tag); }

  //Get the tag associated with an index; this is suboptimal obviously
  TagType tag(const size_t idx) const{ 
    for(auto it=tag_map->begin(); it !=tag_map->end(); it++) if(it->second == idx) return it->first;
    error_exit(std::cout << printType<taggedValueContainer<T> >() << "::tag could not find index " << idx << " in tag map\n");
  }

  //Import/export to std::map
  //Optionally build the tag map using the ordering of the input map
  void mapImport(const std::map<TagType,T> &in, bool generate_tag_map = false){
    if(!generate_tag_map)  assert(tag_map != NULL);
    else{
      static std::list<TagIdxMap> todelete;
      TagIdxMap mp;
      size_t i=0;
      for(auto it=in.begin(); it != in.end(); it++) mp[it->first] = i++;
      todelete.push_back(std::move(mp));
      tag_map = &todelete.back();
    }      
    t.resize(tag_map->size());

    for(auto it = tag_map->begin(); it != tag_map->end(); it++){
      auto in_it = in.find(it->first);
      if(in_it == in.end()) error_exit(std::cout << "taggedValueContainer::import input map does not contain tag " << it->first << std::endl);
      t[it->second] = in_it->second;
    }
  }
  void mapExport(std::map<TagType,T> &out) const{
    out.clear();
    assert(tag_map != NULL);
    for(auto it = tag_map->begin(); it != tag_map->end(); it++){
      out[it->first] = t[it->second];
    }
  }
};

template<typename T, typename TagType>
struct getElem<taggedValueContainer<T,TagType> >{
  static inline T& elem(taggedValueContainer<T,TagType> &v, const int i){ return v(i); }
  static inline const T& elem(const taggedValueContainer<T,TagType> &v, const int i){ return v(i); }    
  inline static typename taggedValueContainer<T,TagType>::TagIdxMap const* common_properties(const taggedValueContainer<T,TagType> &v){ return v.getTagMap(); } //for setting up ET output
};
template<typename T, typename TagType>
std::ostream & operator<<(std::ostream &os, const taggedValueContainer<T,TagType> &r){
  os << r.print(); return os;
}

SARLAC_END_NAMESPACE

#endif
