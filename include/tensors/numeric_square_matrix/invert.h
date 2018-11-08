#ifndef _CPSFIT_NUMERIC_SQUARE_MATRIX_INVERT_H_
#define _CPSFIT_NUMERIC_SQUARE_MATRIX_INVERT_H_

#include<tensors/gsl_svdinverse.h>
#include<tensors/numeric_square_matrix/class.h>
#include<tensors/numeric_square_matrix/sample_view.h>

CPSFIT_START_NAMESPACE

template<typename Numeric, ENABLE_IF_FLOATINGPT(Numeric)>
inline int svd_inverse(NumericSquareMatrix<Numeric> &Ainv, 
		       const NumericSquareMatrix<Numeric> &A,
		       Numeric &condition_number){ 
  return GSL_SVDinvert<NumericSquareMatrix<Numeric>, NumericSquareMatrix<Numeric> >::doit(Ainv,A,condition_number);
}
template<typename Numeric, ENABLE_IF_FLOATINGPT(Numeric)>
inline int svd_inverse(NumericSquareMatrixSampleView<NumericSquareMatrix<Numeric> > &Ainv, 
		       NumericSquareMatrixSampleView<const NumericSquareMatrix<Numeric> > &A,
		       Numeric &condition_number){ 
  return GSL_SVDinvert<NumericSquareMatrixSampleView<NumericSquareMatrix<Numeric> >, 
		       NumericSquareMatrixSampleView<const NumericSquareMatrix<Numeric> > >::doit(Ainv,A,condition_number);
}

//For Matrix<D> where D is a distribution or any type with a 'sample' method
template<typename NumericDist, typename std::enable_if<hasSampleMethod<NumericDist>::value, int>::type = 0>
int svd_inverse(NumericSquareMatrix<NumericDist> &Ainv, 
	        const NumericSquareMatrix<NumericDist> &A,
		NumericDist &condition_number){ 
  const int nsample = A(0,0).size();
  int ret = 0;
  std::vector<int> ret_thr(omp_get_max_threads(),0);
  
  condition_number.resize(nsample);

#pragma omp parallel for
  for(int j=0;j<nsample;j++){
    NumericSquareMatrixSampleView<const NumericSquareMatrix<NumericDist> > A_view(A,j);
    NumericSquareMatrixSampleView<NumericSquareMatrix<NumericDist> > Ainv_view(Ainv,j);

    int r = GSL_SVDinvert<decltype(Ainv_view), decltype(A_view)>::doit(Ainv_view,A_view,condition_number.sample(j));
    ret_thr[omp_get_thread_num()] = ret_thr[omp_get_thread_num()] || r;
  }
  for(int i=0;i<omp_get_max_threads();i++) ret = ret || ret_thr[i];
  
  return ret;
}

template<typename Numeric>
inline int svd_inverse(NumericSquareMatrix<Numeric> &Ainv, 
		       const NumericSquareMatrix<Numeric> &A){
  Numeric c;
  return svd_inverse(Ainv, A, c);
}

CPSFIT_END_NAMESPACE

#endif
