#ifndef _SARLAC_NUMERIC_TENSOR_IO_H_
#define _SARLAC_NUMERIC_TENSOR_IO_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_tensor/class.h>

SARLAC_START_NAMESPACE

template<typename T, int R>
void _NTprinter<T,R>::print(std::ostream &os, const NumericTensor<T,R> &t){
  int e[R];
  for(size_t i=0;i<t.vol;i++){
    t.unmap(e,i);
    os << '(';
    for(int d=0;d<R-1;d++) os << e[d] << ',';
    os << e[R-1] << ") : " << t.data[i] << std::endl;
  }
};

//Matrix specialization
template<typename T>
struct _NTprinter<T,2>{
  static void print(std::ostream &os, const NumericTensor<T,2> &t){
    for(int i=0;i<t.size(0);i++){
      for(int j=0;j<t.size(1)-1;j++){
	os << t({i,j}) << ", ";
      }
      os << t({i,t.size(1)-1}) << "\n";
    }
  }
};
//Vector specialization
template<typename T>
struct _NTprinter<T,1>{
  static void print(std::ostream &os, const NumericTensor<T,1> &t){
    for(int j=0;j<t.size(0)-1;j++){
      os << t({j}) << ", ";
    }
    os << t({t.size(0)-1});
  }
};

template<typename T, int R>
inline std::ostream & operator<<(std::ostream &os, const NumericTensor<T,R> &t){
  _NTprinter<T,R>::print(os,t);
  return os;
}

SARLAC_END_NAMESPACE
#endif
