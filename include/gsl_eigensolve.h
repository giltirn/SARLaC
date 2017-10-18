#ifndef _GSL_EIGENSOLVE_H_
#define _GSL_EIGENSOLVE_H_

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <template_wizardry.h>
#include <numeric_tensors.h>
#include <distribution_iterate.h>

//Requires a floating point *square* matrix with  (i,j) accessor and size() operation
//evecs and evals size should be equal to matrix size, and eigenvector size should too
template<typename VectorOutputType, typename MatrixInputType,
	 typename std::enable_if<
	   std::is_floating_point<typename _get_elem_type<MatrixInputType>::type>::value &&
	   std::is_same<typename _get_elem_type<MatrixInputType>::type, typename _get_vector_elem_type<VectorOutputType>::type>::value
	   , int>::type = 0
	 >
struct _eigensolver{
  static void symmetricMatrixSolve(std::vector<VectorOutputType> &evecs, std::vector<double> &evals, const MatrixInputType &A, bool sort = true){
    typedef typename _get_elem_type<MatrixInputType>::type T;
    const int size = A.size();
    assert(evecs.size() == size);
    assert(evals.size() == size);
    
    gsl_vector *eval = gsl_vector_alloc(size);
    gsl_matrix *evec = gsl_matrix_alloc(size, size);
    gsl_matrix *m  = gsl_matrix_alloc(size,size);
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	gsl_matrix_set(m,i,j,A(i,j));
      }
    }

    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(size);
    gsl_eigen_symmv(m, eval, evec, w);
    gsl_eigen_symmv_free(w);

    if(sort) gsl_eigen_symmv_sort(eval, evec, 
				  GSL_EIGEN_SORT_ABS_ASC);

    for(int i=0;i<size;i++){
      evals[i] = gsl_vector_get(eval,i);
      for(int j=0;j<size;j++)
	evecs[i](j) = gsl_matrix_get(evec,i,j);
    }
    gsl_matrix_free(m);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
  }
};

template<typename T, int is_distribution>
struct _symmetricMatrixEigensolve{};

template<typename T>
struct _symmetricMatrixEigensolve<T,0>{
  inline static void solve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A){
    _eigensolver<NumericVector<T>, NumericSquareMatrix<T> >::symmetricMatrixSolve(evecs,evals,A);
  }
};

template<typename T>
struct _symmetricMatrixEigensolve<T,1>{
  static void solve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A){
    assert(A.size() > 0);
    const int size = A.size();
    const int nsample = A(0,0).size();

    for(int i=0;i<size;i++){
      evals[i].resize(nsample);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(nsample);
    }
    
    typedef typename std::decay< decltype( A(0,0).sample(0) ) >::type type;
    NumericSquareMatrix<type> A_s(size);
    std::vector<NumericVector<type> > evecs_s(size, NumericVector<type>(size));
    std::vector<type> evals_s(size);

    const int nit = iterate<T>::size(A(0,0));
    for(int s=0;s<nit;s++){
      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	
      _eigensolver<NumericVector<type>, NumericSquareMatrix<type> >::symmetricMatrixSolve(evecs_s,evals_s,A_s, false); //don't sort

      for(int i=0;i<size;i++){
	iterate<T>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<T>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
  }
};


template<typename T>
void symmetricMatrixEigensolve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A){
  const int size = A.size();
  evecs.resize(size, NumericVector<T>(size));
  evals.resize(size);
  _symmetricMatrixEigensolve<T, hasSampleMethod<T>::value>::solve(evecs,evals,A);
}


#endif
