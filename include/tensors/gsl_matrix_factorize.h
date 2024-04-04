#pragma once

#include<cassert>
#include<cstring>
#include<gsl/gsl_math.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_linalg.h>

#include<tensors/numeric_vector.h>
#include<config.h>
#include<utils/template_wizardry.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

template<typename MatrixInputType>
struct GSLmatrixFactorize{
  //Compute the Cholesky factorization A=L L^\dagger  for lower-triangular L. Requires symmetric, positive-definite matrix (i.e it is Hermitian)
  template<typename U = MatrixInputType, typename std::enable_if<  is_std_complex<typename _get_elem_type<U>::type>::value, int>::type = 0>
  static void Cholesky(MatrixInputType &L, const MatrixInputType &A){
    const int size = A.size();
    assert(L.size() == size);
    typedef typename _get_elem_type<MatrixInputType>::type complexType;
    gsl_matrix_complex *Ag = gsl_matrix_complex_alloc(size, size);
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	gsl_complex v; 
	GSL_REAL(v) = A(i,j).real();
	GSL_IMAG(v) = A(i,j).imag();
	gsl_matrix_complex_set(Ag,i,j,v);
      }
    }
    int ret = gsl_linalg_complex_cholesky_decomp(Ag);
    if(ret) throw std::runtime_error(std::string("Cholesky (complex) failed with error: ") + std::string(gsl_strerror(ret)));
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	if(i>=j){
	  gsl_complex v = gsl_matrix_complex_get(Ag,i,j);
	  L(i,j) = complexType(GSL_REAL(v),GSL_IMAG(v));
	}else{
	  L(i,j) = complexType(0.);
	}
      }
    }
    gsl_matrix_complex_free(Ag);
  }

  //Compute the Cholesky factorization A=L L^T  for lower-triangular L. Requires symmetric, positive-definite matrix
  template<typename U = MatrixInputType, typename std::enable_if<  std::is_floating_point<typename _get_elem_type<U>::type>::value, int>::type = 0>
  static void Cholesky(MatrixInputType &L, const MatrixInputType &A){
    const int size = A.size();
    assert(L.size() == size);
    gsl_matrix *Ag = gsl_matrix_alloc(size, size);
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	gsl_matrix_set(Ag,i,j,A(i,j));
      }
    }
    int ret = gsl_linalg_cholesky_decomp(Ag);
    if(ret) throw std::runtime_error(std::string("Cholesky (real) failed with error: ") + std::string(gsl_strerror(ret)));
    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	if(i>=j){
	  double v = gsl_matrix_get(Ag,i,j);
	  L(i,j) = v;
	}else{
	  L(i,j) = 0.;
	}
      }
    }
    gsl_matrix_free(Ag);
  }
};


SARLAC_END_NAMESPACE
