#ifndef _SARLAC_PLOT_WRITE_PYTHON_H_
#define _SARLAC_PLOT_WRITE_PYTHON_H_

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>

SARLAC_START_NAMESPACE

//Print out a generic type as a string for use in formulating the Matplotlib kwarg dictionaries
class KWargElem: public OstreamHook{
  std::string val;

public:
  KWargElem(){}
    
  template<typename T>
  KWargElem(const T& v){
    *this = v;
  }

  template<typename T>
  KWargElem &operator=(const T &v){
    std::ostringstream os; os << v;
    val = os.str();
    return *this;
  }

  KWargElem &operator=(const char* v){
    val = "\"" + std::string(v) + "\"";
    return *this;
  }
  
  KWargElem &operator=(const std::string &v){
    val = "\"" + v + "\"";
    return *this;
  }
  KWargElem &operator=(const char v){
    std::ostringstream os; os << "\'" << v << "\'";
    val = os.str();
    return *this;
  }

  KWargElem &operator=(const bool v){
    val = v ? "True" : "False";
    return *this;
  }
  
  inline void write(std::ostream &os) const{
    os << val;
  }
};

//Print a vector as a python list
template<typename T>
struct ListPrint: public OstreamHook{
  const std::vector<T> &lst;
  const std::string elem_enclose;
  ListPrint(const std::vector<T> &_lst, const std::string _elem_enclose = ""): lst(_lst), elem_enclose(_elem_enclose){}
  void write(std::ostream &os) const{
    os << '[';
    for(int i=0;i<lst.size()-1; i++) os << elem_enclose << lst[i] << elem_enclose << ',';      
    os << elem_enclose << lst.back() << elem_enclose << ']';
  }
};

template<typename T>
class PythonTuple{
  T v[2];
public:
  PythonTuple(){}
  PythonTuple(const T&a, const T&b){ v[0]=a; v[1]=b; }

  inline T & operator[](const int i){ return v[i]; }
  inline const T & operator[](const int i) const{ return v[i]; }
};
template<typename T>
std::ostream & operator<<(std::ostream &os, const PythonTuple<T> &v){
  os << "(" << v[0] <<  ", " << v[1] << ")";
  return os;
}


SARLAC_END_NAMESPACE
#endif
