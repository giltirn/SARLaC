#ifndef _SARLAC_NUMERIC_SQUARE_MATRIX_IO_H_
#define _SARLAC_NUMERIC_SQUARE_MATRIX_IO_H_

#include<tensors/numeric_square_matrix/class.h>

SARLAC_START_NAMESPACE

#ifdef HAVE_HDF5
template<typename D>
inline void write(HDF5writer &writer, const NumericSquareMatrix<D> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, NumericSquareMatrix<D> &d, const std::string &tag){ d.read(reader,tag); }
#endif

template<typename Numeric, typename StreamType, typename std::enable_if< isStreamType<StreamType>::value, int>::type = 0> 
StreamType & operator<<(StreamType & stream, const NumericSquareMatrix<Numeric> &mat){
  for(int i=0;i<mat.size();i++){
    for(int j=0;j<mat.size();j++){
      stream << mat(i,j) << " ";
    }
    stream << '\n';
  }
  return stream;
}

SARLAC_END_NAMESPACE
#endif
