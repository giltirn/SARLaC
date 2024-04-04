#pragma once

#include<tensors/numeric_rect_matrix/class.h>

SARLAC_START_NAMESPACE

#ifdef HAVE_HDF5
template<typename D>
inline void write(HDF5writer &writer, const NumericRectMatrix<D> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, NumericRectMatrix<D> &d, const std::string &tag){ d.read(reader,tag); }
#endif

template<typename Numeric, typename StreamType, typename std::enable_if< isStreamType<StreamType>::value, int>::type = 0> 
StreamType & operator<<(StreamType & stream, const NumericRectMatrix<Numeric> &mat){
  for(int i=0;i<mat.rows();i++){
    for(int j=0;j<mat.cols();j++){
      stream << mat(i,j) << " ";
    }
    stream << '\n';
  }
  return stream;
}

SARLAC_END_NAMESPACE
