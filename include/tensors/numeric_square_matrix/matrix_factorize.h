#pragma once

#include<mutex>
#include<distribution/distribution_iterate.h>
#include<tensors/gsl_matrix_factorize.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix/class.h>
#include<utils/template_wizardry/complexify.h>
#include<utils/template_wizardry/number.h>

SARLAC_START_NAMESPACE

//Perform a Cholesky decomposition A = L L^\dagger for lower-triangular L  for symmetric, positive-definite A. Returns L
template<typename T, typename std::enable_if< is_floating_point_or_complex<T>::value, int>::type = 0>
NumericSquareMatrix<T> Cholesky(const NumericSquareMatrix<T> &A){
  NumericSquareMatrix<T> L(A.size());
  GSLmatrixFactorize<NumericSquareMatrix<T> >::Cholesky(L,A);
  return L;
}

SARLAC_END_NAMESPACE
