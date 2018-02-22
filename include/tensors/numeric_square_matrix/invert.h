#ifndef _CPSFIT_NUMERIC_SQUARE_MATRIX_INVERT_H_
#define _CPSFIT_NUMERIC_SQUARE_MATRIX_INVERT_H_

#include<tensors/gsl_svdinverse.h>
#include<tensors/numeric_square_matrix/class.h>
#include<tensors/numeric_square_matrix/sample_view.h>

CPSFIT_START_NAMESPACE

template<typename NumericSquareMatrixOutputType, typename NumericSquareMatrixInputType, ENABLE_IF_ELEM_TYPE_FLOATINGPT(NumericSquareMatrixInputType)>
int svd_inverse(NumericSquareMatrixOutputType &Ainv, 
	        const NumericSquareMatrixInputType &A){ 
  double c;
  return GSL_SVDinvert<NumericSquareMatrixOutputType, NumericSquareMatrixInputType>::doit(Ainv,A,c);
}


//For Matrix<D> where D is a distribution or any type with a 'sample' method
template<typename NumericSquareMatrixType, typename std::enable_if<hasSampleMethod<typename _get_elem_type<NumericSquareMatrixType>::type>::value, int>::type = 0>
int svd_inverse(NumericSquareMatrixType &Ainv, 
	        const NumericSquareMatrixType &A){ 
  const int nsample = A(0,0).size();
  int ret = 0;
  std::vector<int> ret_thr(omp_get_max_threads(),0);
  
#pragma omp parallel for
  for(int j=0;j<nsample;j++){
    NumericSquareMatrixSampleView<const NumericSquareMatrixType> A_view(A,j);
    NumericSquareMatrixSampleView<NumericSquareMatrixType> Ainv_view(Ainv,j);
    int r = svd_inverse(Ainv_view,A_view);
    ret_thr[omp_get_thread_num()] = ret_thr[omp_get_thread_num()] || r;
  }
  for(int i=0;i<omp_get_max_threads();i++) ret = ret || ret_thr[i];
  
  return ret;
}


CPSFIT_END_NAMESPACE

#endif
