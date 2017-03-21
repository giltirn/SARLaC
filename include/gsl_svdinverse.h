#ifndef _GSL_SVDINVERSE_H_
#define _GSL_SVDINVERSE_H_

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <template_wizardry.h>

template<typename MatrixType>
struct _get_elem_type{
  typedef decltype( ( (MatrixType*)(nullptr) )->operator()(0,0) ) RefType;
  typedef typename std::remove_reference<RefType>::type BaseType;
  typedef typename std::remove_const<BaseType>::type type;
};

//Requires a floating point *square* matrix with  (i,j) accessor and size() operation
template<typename MatrixOutputType, typename MatrixInputType,
	 typename std::enable_if<
	   std::is_floating_point<typename _get_elem_type<MatrixInputType>::type>::value
	   , int>::type = 0
	 >
struct _svd_inverse{
  static int doit(MatrixOutputType &Ainv, const MatrixInputType &A, double &condition_number){
    typedef typename _get_elem_type<MatrixInputType>::type T;
    
    int size = A.size();
    assert(Ainv.size() == size);
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
  int size() const{ return M.size(); }

  T_const& operator()(const int i, const int j) const{ return M[i][j]; }
  
  template<typename U = T>
  typename std::enable_if< !std::is_const<U>::value, T_unconst& >::type operator()(const int i, const int j){ return M[i][j]; } //unconst accessor only for unconst T
};

template<typename T>  
int svd_inverse(std::vector< std::vector<T> > &Ainv, 
	        const std::vector< std::vector<T> > &A){
  vectorVectorMatrixView<T> Ainv_view(Ainv);
  vectorVectorMatrixView<const T> A_view(A);  
  double c;
  return _svd_inverse<vectorVectorMatrixView<T>, vectorVectorMatrixView<const T> >::doit(Ainv_view,A_view,c);
}
template<typename T>  
int svd_inverse(std::vector< std::vector<T> > &Ainv, 
	        const std::vector< std::vector<T> > &A,
		double &condition_number){
  vectorVectorMatrixView<T> Ainv_view(Ainv);
  vectorVectorMatrixView<const T> A_view(A);  
  return _svd_inverse<vectorVectorMatrixView<T>, vectorVectorMatrixView<const T> >::doit(Ainv_view,A_view,condition_number);
}




#endif
