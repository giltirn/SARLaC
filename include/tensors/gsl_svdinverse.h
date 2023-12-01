#ifndef _GSL_SVDINVERSE_H_
#define _GSL_SVDINVERSE_H_

//Real square matrix inverse using GSL SVD algorithm

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>

SARLAC_START_NAMESPACE

//Requires a *real* floating point *square* matrix with  (i,j) accessor and size() operation
template<typename MatrixOutputType, typename MatrixInputType,
	 ENABLE_IF_ELEM_TYPE_FLOATINGPT(MatrixInputType)
	 >
struct GSL_SVDinvert{
  static int doit(MatrixOutputType &Ainv, const MatrixInputType &A, double &condition_number){
    typedef typename _get_elem_type<MatrixInputType>::type T;
    
    int size = A.size();
    assert(Ainv.size() == size);
    if(size == 1){
      Ainv(0,0) = 1.0/A(0,0);
      return 0;
    }
    int ret;

    std::vector< std::vector<T> > U(size, std::vector<T>(size));
    std::vector< std::vector<T> > V(size, std::vector<T>(size));
    std::vector<T> Diag(size);

    gsl_matrix *m  = gsl_matrix_alloc(size,size);
    gsl_matrix *Q  = gsl_matrix_alloc(size,size);
    gsl_vector *S  = gsl_vector_alloc(size);
    gsl_vector *WORK  = gsl_vector_alloc(size);

    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	gsl_matrix_set(m,i,j,A(i,j));
      }
    }
    //ret = gsl_linalg_SV_decomp_jacobi(m,Q,S);
    ret = gsl_linalg_SV_decomp(m,Q,S,WORK);

    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	U[i][j] = gsl_matrix_get(m,i,j);
	V[i][j] = gsl_matrix_get(Q,i,j);
	Ainv(i,j) = 0.0;
      }
    }
    for(int i=0;i<size;i++){
      Diag[i] = 1.0/gsl_vector_get(S,i); 
    }     
    int smallest_sing_idx=size-1;

    //When singular value sigma_i << sigma_max, numerical errors can lead to large errors in the inverted matrix.
    //It is therefore common practise to discard these singular values.

  // if(FitGlobals::svd_singularity_threshold >0.0){
  //   for(int i=size-1;i>=0;i--){
  //     if(gsl_vector_get(S,i)/gsl_vector_get(S,0)<=FitGlobals::svd_singularity_threshold){ 
  // 	Diag[i] = 0.0; --smallest_sing_idx; 
  // 	printf("SVD sing threshold discard %d\n",i);
  //     }
  //     else break; //singular values are automatically ordered by gsl routine, ;argest first
  //   }   
  // }

    for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
	for(int k=0;k<size;k++){
	  Ainv(i,k) += V[i][j]*Diag[j]*U[k][j];  
	}}}

    condition_number = gsl_vector_get(S,0)/gsl_vector_get(S,smallest_sing_idx);
    
    
    gsl_matrix_free(m);
    gsl_matrix_free(Q);
    gsl_vector_free(S);
    gsl_vector_free(WORK);
    
    if ( ret ) { 
      printf("Warning GSL inversion failed %d\n",ret);
    }
    return ret;    
  }

};

template<typename T>
class vectorVectorMatrixView{
  typedef typename std::remove_const<T>::type T_unconst;
  typedef typename std::add_const<T_unconst>::type T_const;
  typedef typename std::vector<std::vector<T_unconst> > MtypeBase;
  typedef typename add_const_if<MtypeBase, T>::type Mtype;   //const MtypeBase or MtypeBase if T const or unconst respectively
  
  Mtype &M;
public:
  vectorVectorMatrixView(Mtype &_M): M(_M){ assert(M[0].size() == M.size()); }
  inline int size() const{ return M.size(); }

  inline T_const& operator()(const int i, const int j) const{ return M[i][j]; }
  
  template<typename U = T>
  inline typename std::enable_if< !std::is_const<U>::value, T_unconst& >::type operator()(const int i, const int j){ return M[i][j]; } //unconst accessor only for unconst T
};

template<typename T>  
int svd_inverse(std::vector< std::vector<T> > &Ainv, 
	        const std::vector< std::vector<T> > &A){
  vectorVectorMatrixView<T> Ainv_view(Ainv);
  vectorVectorMatrixView<const T> A_view(A);  
  double c;
  return GSL_SVDinvert<vectorVectorMatrixView<T>, vectorVectorMatrixView<const T> >::doit(Ainv_view,A_view,c);
}
template<typename T>  
int svd_inverse(std::vector< std::vector<T> > &Ainv, 
	        const std::vector< std::vector<T> > &A,
		double &condition_number){
  vectorVectorMatrixView<T> Ainv_view(Ainv);
  vectorVectorMatrixView<const T> A_view(A);  
  return GSL_SVDinvert<vectorVectorMatrixView<T>, vectorVectorMatrixView<const T> >::doit(Ainv_view,A_view,condition_number);
}


//Compute the Moore-Penrose pseudo-inverse (based on https://gist.github.com/turingbirds/5e99656e08dbe1324c99)
//Requires a floating point matrix with  (i,j) accessor,  int rows() and int cols()
//rcond is the threshold singular value below which it is considered to be exactly zero
template<typename MatrixOutputType, typename MatrixInputType,
	 ENABLE_IF_ELEM_TYPE_FLOATINGPT(MatrixInputType)
	 >
struct GSL_MoorePenrosePseudoInverse{
  static void doit(MatrixOutputType &Ainv_, const MatrixInputType &A_, const double rcond = 1e-15){
    typedef typename _get_elem_type<MatrixInputType>::type T;

    size_t n = A_.rows();
    size_t m = A_.cols();

    gsl_matrix *A_base = gsl_matrix_alloc(n, m);
    
    gsl_matrix *A = A_base;
    for(size_t i=0;i<n;i++)
      for(size_t j=0;j<m;j++)
	gsl_matrix_set(A, i, j, A_(i,j));
    
    gsl_matrix *_tmp_mat = NULL;
    double x, cutoff;
    bool was_swapped = false;

    if (m > n) {
      /* libgsl SVD can only handle the case m <= n - transpose matrix */
      was_swapped = true;
      _tmp_mat = gsl_matrix_alloc(m, n);
      gsl_matrix_transpose_memcpy(_tmp_mat, A);
      A = _tmp_mat;
      size_t i = m;
      m = n;
      n = i;
    }

    /* do SVD */
    gsl_matrix *V = gsl_matrix_alloc(m, m);
    gsl_vector *u = gsl_vector_alloc(m);
    gsl_vector *_tmp_vec = gsl_vector_alloc(m);
    gsl_linalg_SV_decomp(A, V, u, _tmp_vec);
    gsl_vector_free(_tmp_vec);

    /* compute \Sigma^{-1} */
    gsl_matrix *Sigma_pinv = gsl_matrix_alloc(m, n);
    gsl_matrix_set_zero(Sigma_pinv);
    cutoff = rcond * gsl_vector_max(u);
    
    for (size_t i = 0; i < m; ++i) {
      x = gsl_vector_get(u, i) > cutoff ? 1. / gsl_vector_get(u, i) : 0.;
      gsl_matrix_set(Sigma_pinv, i, i, x);
    }
    
    /* libgsl SVD yields "thin" SVD - pad to full matrix by adding zeros */
    gsl_matrix *U = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(U);
    
    for(size_t i = 0; i < n; ++i)
      for(size_t j = 0; j < m; ++j)
	gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));

    if(_tmp_mat != NULL)
      gsl_matrix_free(_tmp_mat);

    /* two dot products to obtain pseudoinverse */
    _tmp_mat = gsl_matrix_alloc(m, n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Sigma_pinv, 0., _tmp_mat);
    
    gsl_matrix *A_pinv;
    if (was_swapped) {
      A_pinv = gsl_matrix_alloc(n, m);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., U, _tmp_mat, 0., A_pinv);

      for(int i=0;i<n;i++)
	for(int j=0;j<m;j++)
	  gsl_matrix_set(A, i, j, A_(i,j));
    }
    else {
      A_pinv = gsl_matrix_alloc(m, n);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., _tmp_mat, U, 0., A_pinv);

      for(int i=0;i<m;i++)
	for(int j=0;j<n;j++)
	  gsl_matrix_set(A, i, j, A_(i,j));
    }
    
    gsl_matrix_free(A_base);
    gsl_matrix_free(_tmp_mat);
    gsl_matrix_free(U);
    gsl_matrix_free(Sigma_pinv);
    gsl_vector_free(u);
    gsl_matrix_free(V);
    gsl_matrix_free(A_pinv);
  }
};








SARLAC_END_NAMESPACE

#endif
