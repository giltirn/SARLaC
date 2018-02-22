#ifndef _CPSFIT_NUMERIC_VECTOR_IO_H_
#define _CPSFIT_NUMERIC_VECTOR_IO_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector/class.h>

CPSFIT_START_NAMESPACE

#ifdef HAVE_HDF5
template<typename D>
inline void write(HDF5writer &writer, const NumericVector<D> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, NumericVector<D> &d, const std::string &tag){ d.read(reader,tag); }
#endif


template<typename Numeric> 
std::ostream & operator<<(std::ostream & stream, const NumericVector<Numeric> &vec){
  stream << "(";
  for(int i=0;i<vec.size();i++)
    stream << vec[i] << (i != vec.size()-1 ? " " : ")");
  return stream;
}

CPSFIT_END_NAMESPACE

#endif
