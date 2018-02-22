#ifndef _GSL_EIGENSOLVE_H_
#define _GSL_EIGENSOLVE_H_

//A wrapper around GSL eigensolver for symmetric matrices
#include<cassert>
#include<gsl/gsl_math.h>
#include<gsl/gsl_eigen.h>

#include<config.h>
#include<utils/template_wizardry.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//Requires a real floating point *square* matrix with  (i,j) accessor and size() operation
//evecs and evals size should be equal to matrix size, and eigenvector size should too
template<typename VectorOutputType, typename MatrixInputType,
	 typename std::enable_if<
	   std::is_floating_point<typename _get_elem_type<MatrixInputType>::type>::value &&
	   std::is_same<typename _get_elem_type<MatrixInputType>::type, typename _get_vector_elem_type<VectorOutputType>::type>::value
	   , int>::type = 0
	 >
struct GSLsymmEigenSolver{
  static std::vector<double> symmetricMatrixSolve(std::vector<VectorOutputType> &evecs, std::vector<double> &evals, const MatrixInputType &A, bool sort = true){
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
				  GSL_EIGEN_SORT_VAL_DESC);

    for(int i=0;i<size;i++){
      evals[i] = gsl_vector_get(eval,i);
      for(int j=0;j<size;j++)
	evecs[i](j) = gsl_matrix_get(evec,j,i); //eigenvectors stored in columns, i.e. elements have fixed column index and changing row index
    }
    gsl_matrix_free(m);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);

    std::vector<double> residuals(size);
    for(int i=0;i<size;i++){
      VectorOutputType Mv_m_lv = A * evecs[i] - evals[i] * evecs[i];
      residuals[i] = sqrt(mod2(Mv_m_lv));
    }
    return residuals;
  }
};

CPSFIT_END_NAMESPACE

#endif
