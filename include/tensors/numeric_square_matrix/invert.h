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
  const int nsample = iterate<NumericDist>::size(A(0,0));

  const static std::string fname = "svd_inverse(NumericSquareMatrix<NumericDist> &Ainv, const NumericSquareMatrix<NumericDist> &A, NumericDist &condition_number)";

  if(Ainv.size() != A.size()) error_exit(std::cout << fname << " inverse matrix must have same size as input matrix\n");

  for(int i=0;i<A.size();i++){
    for(int j=0;j<A.size();j++){
      if(iterate<NumericDist>::size(A(i,j)) != nsample) 
	error_exit(std::cout << fname << " sample number discrepancy for A(" << i << ","<< j << ") got " << iterate<NumericDist>::size(A(i,j)) << " expect " << nsample << "\n");
      if(iterate<NumericDist>::size(Ainv(i,j)) != nsample) 
	error_exit(std::cout << fname << " sample number discrepancy for Ainv(" << i << ","<< j << ") got " << iterate<NumericDist>::size(Ainv(i,j)) << " expect " << nsample << "\n");
    }
  }

  int ret = 0;
  std::vector<int> ret_thr(omp_get_max_threads(),0);
  
  condition_number = A(0,0);

#pragma omp parallel for
  for(int j=0;j<nsample;j++){
    NumericSquareMatrixSampleView<const NumericSquareMatrix<NumericDist> > A_view(A,j);
    NumericSquareMatrixSampleView<NumericSquareMatrix<NumericDist> > Ainv_view(Ainv,j);

    int r = GSL_SVDinvert<decltype(Ainv_view), decltype(A_view)>::doit(Ainv_view,A_view, iterate<NumericDist>::at(j, condition_number));
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
