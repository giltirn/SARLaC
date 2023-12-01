#ifndef _SET_OPERATIONS_H_
#define _SET_OPERATIONS_H_

#include<set>
#include<vector>

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE


//both sets together
template<typename T>
inline std::set<T> setUnion(const std::set<T> &a, const std::set<T> &b){
  std::set<T> out;
  for(auto it=a.begin(); it!=a.end(); ++it) out.insert(*it);
  for(auto it=b.begin(); it!=b.end(); ++it) out.insert(*it);
  return out;
}

//members of both a and b
template<typename T>
inline std::set<T> setIntersection(const std::set<T> &a, const std::set<T> &b){
  std::set<T> out;
  for(auto it=a.begin(); it!=a.end(); ++it) if(b.count(*it)) out.insert(*it);
  return out;
}

//members of a that are not members of b
template<typename T>
inline std::set<T> setComplement(const std::set<T> &a, const std::set<T> &b){
  std::set<T> out;
  for(auto it=a.begin(); it!=a.end(); ++it) if(!b.count(*it)) out.insert(*it);
  return out;
}

template<typename T>
inline std::vector<T> setToVector(const std::set<T> &a){
  std::vector<T> out(a.size()); 
  size_t i = 0; 
  for(auto it=a.begin(); it!=a.end(); ++it) out[i++] = *a;
  return out;
}


SARLAC_END_NAMESPACE

#endif
