#ifndef _CPSFIT_NUMERIC_SQUARE_MATRIX_MATH_H_
#define _CPSFIT_NUMERIC_SQUARE_MATRIX_MATH_H_

#include<tensors/numeric_square_matrix/class.h>

CPSFIT_START_NAMESPACE

//Modulus operator for matrices
template<typename T>
T mod2(const NumericSquareMatrix<T> &m){
  assert(m.size() > 0);
  T out(m(0,0));
  zeroit(out);
  for(int i=0;i<m.size();i++)
    for(int j=0;j<m.size();j++)
      out = out + m(i,j)*m(i,j); //only correct for real matrices!
  return out;
}

CPSFIT_END_NAMESPACE
#endif
