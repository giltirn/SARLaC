#ifndef _CPSFIT_UTILS_STRING_H_
#define _CPSFIT_UTILS_STRING_H_

//String manipulation
#include<iostream>
#include<map>
#include<cassert>
#include<sstream>
#include<cstdarg>
#include<string>
#include<algorithm>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils/error.h>

CPSFIT_START_NAMESPACE

//A bi-directional mapping between a string tag and an integer
class stringTagMap{
  std::map<std::string,int> mp;
  std::map<int,std::string> ump;
public:
  typedef std::string tagType; 
  inline int map(const tagType &tag) const{
    auto it = mp.find(tag);
    assert(it != mp.end());
    return it->second;
  }
  inline tagType unmap(const int idx) const{
    auto it = ump.find(idx);
    assert(it != ump.end());
    return it->second;
  }
  void add(const std::string &tag, const int idx){
    mp[tag] = idx;
    ump[idx] = tag;
  }
  inline int size() const{ return mp.size(); }
};


//Substitute substring '%d' with 'idx'
inline std::string subsIdx(const std::string fmt, const int idx){
  std::string::size_type off = fmt.find("%d");
  if(off == std::string::npos){
    std::cout << "Could not find substring \"%d\" in format string " << fmt << std::endl;
    std::cout.flush();
    exit(-1);
  }
  std::ostringstream os; os << idx;
  std::string out(fmt);
  out.replace(off,2,os.str());
  return out;
}

//Convert back and forth between string and a general type
template<typename T>
inline T strToAny(const std::string &str){
  T out;
  std::stringstream os(str); os >> out;
  return out;
}
template<typename T>
inline std::string anyToStr(const T &p){
  std::ostringstream os; os << p; return os.str();
}

//C-style string formatting but without the nasty mem buffer concerns
inline std::string stringize(const char* format, ...){
  int n; //not counting null character
  {
    char buf[1024];
    va_list argptr;
    va_start(argptr, format);    
    n = vsnprintf(buf, 1024, format, argptr);
    va_end(argptr);
    if(n < 1024) return std::string(buf);
  }
  
  char buf[n+1];
  va_list argptr;
  va_start(argptr, format);    
  int n2 = vsnprintf(buf, 1024, format, argptr);
  va_end(argptr);
  assert(n2 <= n);
  return std::string(buf);
}


//This class breaks up a string according to a list of tagged substrings that can then be subsequently replaced with new data. 
//This is useful for example in easily handling user-input file format strings

struct subStringSpecify{
  std::string substr;
  bool optional;

  subStringSpecify(std::string substr, bool optional = false): substr(substr),optional(optional){}
};

class subStringReplace{
  std::vector<std::string> chunked;
  std::vector<int> substr_chunk_map; //optional substrings index -1
  std::vector<int> chunk_substr_map;
public:
  void debugPrint() const{
    std::cout << "Chunks:\n";
    for(int i=0;i<chunked.size();i++) std::cout << chunked[i] << std::endl;
    std::cout << "Substring chunk indices:";
    for(int i=0;i<substr_chunk_map.size();i++) std::cout << " " << substr_chunk_map[i];
    std::cout << std::endl;
  }

  bool foundSubstring(const int i) const{ return substr_chunk_map[i] != -1; }

  //Break up 'the_string' around and including the input "tagged_substrings"
  void chunkString(const std::string &the_string, const std::vector<subStringSpecify> &tagged_substrings){
    substr_chunk_map.resize(tagged_substrings.size());
    for(int i=0;i<substr_chunk_map.size();i++) substr_chunk_map[i] = -1;

    std::vector<std::pair<int, size_t> > offsets;
    
    for(int i=0;i<tagged_substrings.size();i++){
      size_t off = the_string.find(tagged_substrings[i].substr);
      if(off == std::string::npos && !tagged_substrings[i].optional) 
	error_exit(std::cout << "Could not find substring \"" << tagged_substrings[i].substr << "\" in " << the_string << std::endl);
      
      if(off != std::string::npos) offsets.push_back(std::pair<int,size_t>(i,off));
    }

    std::sort(offsets.begin(), offsets.end(), [&](const std::pair<int,size_t> &a, const std::pair<int,size_t> &b){ return a.second < b.second; });

    assert(offsets.size() > 0);

    int s = 0;
    size_t pos = 0;

    for(int s=0;s<offsets.size();s++){
      int ssidx = offsets[s].first;
      size_t ssoff = offsets[s].second;

      if(pos != ssoff){
	chunked.push_back(the_string.substr(pos, ssoff-pos));
	pos = ssoff;
      }
    
      chunked.push_back(the_string.substr(pos, tagged_substrings[ssidx].substr.size()));
      pos += tagged_substrings[s].substr.size();
      substr_chunk_map[ssidx] = chunked.size()-1;
    }

    if(pos < the_string.size()-1)
      chunked.push_back(the_string.substr(pos, std::string::npos));

    chunk_substr_map.resize(chunked.size());
    for(int c=0;c<chunked.size();c++) chunk_substr_map[c] = -1;
    
    for(int s=0;s<substr_chunk_map.size();s++) if(substr_chunk_map[s] != -1) chunk_substr_map[substr_chunk_map[s]] = s;
  }

  //'with' should be a vector of strings, one for each substring (including optional)
  void replace(std::ostream &os, const std::vector<std::string> &with) const{
    if(with.size() != substr_chunk_map.size()) error_exit(std::cout << "subStringReplace::replace wrong number of replacement strings provided!\n"); 
    for(int c=0;c<chunked.size();c++){
      int s = chunk_substr_map[c];
      if(s == -1) os << chunked[c];
      else os << with[s];
    }
  }

  subStringReplace() = default;
  
  subStringReplace(const std::string &the_string, const std::vector<subStringSpecify> &tagged_substrings){
    chunkString(the_string, tagged_substrings);
  }
 
};


//Return input string if types match, otherwise return an empty string
template<typename T, typename U>
struct printOnlyIfType{
  inline static std::string str(const std::string &str){ return ""; }
};
template<typename T>
struct printOnlyIfType<T,T>{
  inline static std::string str(const std::string &str){ return str; }
};


CPSFIT_END_NAMESPACE
#endif
