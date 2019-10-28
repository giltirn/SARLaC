#ifndef _GSL_EIGENSOLVE_H_
#define _GSL_EIGENSOLVE_H_

//A wrapper around GSL eigensolver for symmetric matrices
#include<cassert>
#include<cstring>
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


  //GEVP for real, symmetric matrices A, B   solve  A v = lambda B v
  //If there is an error it will throw (i.e. catch if you want to be able to deal with failures)
  static std::vector<double> symmetricGEVPsolve(std::vector<VectorOutputType> &evecs, std::vector<double> &evals, const MatrixInputType &A, const MatrixInputType &B, bool sort = true){
    typedef typename _get_elem_type<MatrixInputType>::type T;
    const int size = A.size();
    assert(B.size() == size);
    assert(evecs.size() == size);
    assert(evals.size() == size);
    
    gsl_vector *eval = gsl_vector_alloc(size);
    gsl_matrix *evec = gsl_matrix_alloc(size, size);
    gsl_matrix *Agsl  = gsl_matrix_alloc(size,size);
    gsl_matrix *Bgsl  = gsl_matrix_alloc(size,size);
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	gsl_matrix_set(Agsl,i,j,A(i,j));
	gsl_matrix_set(Bgsl,i,j,B(i,j));
      }
    }

    gsl_error_handler_t *errh = gsl_set_error_handler_off();

    gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc(size);
    int err = gsl_eigen_gensymmv(Agsl, Bgsl, eval, evec, w);
    gsl_eigen_gensymmv_free(w);

    if(err !=0){
      gsl_set_error_handler(errh); 
      gsl_matrix_free(Agsl);
      gsl_matrix_free(Bgsl);
      gsl_vector_free(eval);
      gsl_matrix_free(evec);
      throw std::runtime_error(gsl_strerror(err));
    }else{
      gsl_set_error_handler (errh); 
      
      if(sort) gsl_eigen_symmv_sort(eval, evec, 
				    GSL_EIGEN_SORT_VAL_DESC);
      
      for(int i=0;i<size;i++){
	evals[i] = gsl_vector_get(eval,i);
	for(int j=0;j<size;j++)
	  evecs[i](j) = gsl_matrix_get(evec,j,i); //eigenvectors stored in columns, i.e. elements have fixed column index and changing row index
      }
      gsl_matrix_free(Agsl);
      gsl_matrix_free(Bgsl);
      gsl_vector_free(eval);
      gsl_matrix_free(evec);
      
      std::vector<double> residuals(size);
      for(int i=0;i<size;i++){
	VectorOutputType Av_m_lBv = A * evecs[i] - evals[i] * B * evecs[i];
	residuals[i] = sqrt(mod2(Av_m_lBv));
      }
      return residuals;
    }
  }

};







//Requires a real floating point *square* matrix with  (i,j) accessor and size() operation
//evecs and evals size should be equal to matrix size, and eigenvector size should too
//Vector output element type should be *complex*
template<typename VectorOutputType, typename MatrixInputType,
	 typename std::enable_if<
	   std::is_floating_point<typename _get_elem_type<MatrixInputType>::type>::value &&
	   std::is_same<std::complex<typename _get_elem_type<MatrixInputType>::type>, typename _get_vector_elem_type<VectorOutputType>::type >::value
	   , int>::type = 0
	 >
struct GSLnonSymmEigenSolver{
  static std::vector<double> nonSymmetricMatrixSolve(std::vector<VectorOutputType> &evecs, std::vector<std::complex<double> > &evals, 
						     const MatrixInputType &A, bool sort = true){
    typedef typename _get_elem_type<MatrixInputType>::type T;
    const int size = A.size();
    assert(evecs.size() == size);
    assert(evals.size() == size);
    
    gsl_vector_complex *eval = gsl_vector_complex_alloc(size);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(size, size);
    gsl_matrix *m  = gsl_matrix_alloc(size,size);
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	gsl_matrix_set(m,i,j,A(i,j));
      }
    }

    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(size);
    gsl_eigen_nonsymmv_params(0,w);

    gsl_eigen_nonsymmv(m, eval, evec, w);
    gsl_eigen_nonsymmv_free(w);

    if(sort) gsl_eigen_nonsymmv_sort(eval, evec, 
				     GSL_EIGEN_SORT_ABS_DESC); //absolute value, descending

    for(int i=0;i<size;i++){
      gsl_complex v = gsl_vector_complex_get(eval,i);
      memcpy(&evals[i], &v, 2*sizeof(T));
      for(int j=0;j<size;j++){
	v = gsl_matrix_complex_get(evec,j,i); //eigenvectors stored in columns, i.e. elements have fixed column index and changing row index
	memcpy(&evecs[i](j), &v, 2*sizeof(T));
      }
    }
    gsl_matrix_free(m);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);

    std::vector<double> residuals(size);
    for(int i=0;i<size;i++){
      VectorOutputType Mv_m_lv = A * evecs[i] - evals[i] * evecs[i];
      double r = 0.;
      for(int j=0;j<Mv_m_lv.size();j++) r += std::abs(Mv_m_lv(j));
      residuals[i] = sqrt(r);
    }
    return residuals;
  }
};



CPSFIT_END_NAMESPACE

#endif
