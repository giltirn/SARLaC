#ifndef _SARLAC_NUMERIC_TENSOR_MATH_H_
#define _SARLAC_NUMERIC_TENSOR_MATH_H_

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_tensor/class.h>
#include<tensors/numeric_square_matrix/invert.h>

SARLAC_START_NAMESPACE

//Contract two tensors over a single index. Preserves ordering of remaining indices
//eg A_ijkl B_mnio -> C_jklmno

template<typename T, int R1, int R2>
NumericTensor<T,R1+R2-2> contract(const NumericTensor<T,R1> &A, const NumericTensor<T,R2> &B, const int rankA, const int rankB){
  if(A.size(rankA) != B.size(rankB)) error_exit(std::cout << "NumericTensor contract(...) contracting over ranks with different sizes!\n");
  const int conRankSize = A.size(rankA);
  int dimC[R1+R2-2];
  int c=0;
  for(int i=0;i<R1;i++)
    if(i!=rankA) dimC[c++] = A.size(i);
  for(int i=0;i<R2;i++)
    if(i!=rankB) dimC[c++] = B.size(i);
  
  return NumericTensor<T,R1+R2-2>(dimC,
				  [&](const int* cC){
				    int cA[R1], cB[R2];
				    int c=0;
				    for(int i=0;i<R1;i++)
				      if(i!=rankA) cA[i]=cC[c++];
				    for(int i=0;i<R2;i++)
				      if(i!=rankB) cB[i]=cC[c++];

				    T out; zeroit(out);
				    for(int i=0;i<conRankSize;i++){
				      cA[rankA] = cB[rankB] = i;				      
				      out = out + A(cA)*B(cB);
				    }
				    return out;
				  });
}


//Functor for tensor reduction by average
template<typename T, int Rank>
struct averageDimensionFunctor{
  const int dim;
  std::vector<int> const* use;
  
  averageDimensionFunctor(const int _dim, std::vector<int> const*_use = NULL): dim(_dim), use(_use){}
  
  void operator()(T &o, int const *coord, const NumericTensor<T,Rank> &from) const{
    int full_coord[Rank];
    int i=0; for(int ii=0;ii<Rank;ii++) if(ii!=dim) full_coord[ii] = coord[i++];    
    zeroit(o);
    if(use != NULL){
      assert(use->size()> 0);
      full_coord[dim] = use->at(0);
      o = from(full_coord);      
      for(int i=1;i<use->size();i++){
	full_coord[dim] = use->at(i);
	o = o + from(full_coord);
      }
      o = o/double(use->size());
    }else{
      assert(from.size(dim)>0);
      full_coord[dim] = 0;
      o = from(full_coord);
      for(int i=1;i<from.size(dim);i++){
	full_coord[dim] = i;
	o = o + from(full_coord);
      }
      o = o/double(from.size(dim));
    }
  }
};
    
template<typename Resampled, typename Raw>
struct resampleFunctor{
  inline Resampled operator()(int const* coord,const Raw &from) const{
    Resampled o(from.size());
    o.resample(from);
    return o;
  }
};

template<typename T>
int svd_inverse(NumericTensor<T,2> &Ainv,
		const NumericTensor<T,2> &A,
		T &condition_number){
  if(A.size(0) != A.size(1)) error_exit(std::cout << "svd_inverse of rank-2 numeric tensor requires a square matrix!\n");
  int N = A.size(0);
  Ainv.resize({N,N});

  NumericSquareMatrix<T> Ainvs(N,A({0,0}));
  NumericSquareMatrix<T> As(N, [&](const int i, const int j){ return A({i,j}); });
  int out = svd_inverse(Ainvs,As,condition_number);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      Ainv({i,j}) = std::move(Ainvs(i,j));
  return out;
}
template<typename T>
inline int svd_inverse(NumericTensor<T,2> &Ainv,
		       const NumericTensor<T,2> &A){
  T c;
  return svd_inverse(Ainv,A,c);
}

template<typename T, ENABLE_IF_FLOATINGPT(T)>
void MoorePenrosePseudoInverse(NumericTensor<T,2> &Ainv,
			       const NumericTensor<T,2> &A,
			       const double rcond = 1e-15){
  struct view{
    NumericTensor<T,2> &M;
    inline size_t rows() const{ return M.size(0); }
    inline size_t cols() const{ return M.size(1); }
    T &operator()(const size_t i, const size_t j){ return M({i,j}); }
    view(NumericTensor<T,2> &M): M(M){}
  };
  struct const_view{
    const NumericTensor<T,2> &M;
    inline size_t rows() const{ return M.size(0); }
    inline size_t cols() const{ return M.size(1); }
    const T &operator()(const size_t i, const size_t j) const{ return M({i,j}); }
    const_view(const NumericTensor<T,2> &M): M(M){}
  };

  Ainv.resize(A.size(1), A.size(0));
  
  const_view in(A);
  view out(Ainv);

  GSL_MoorePenrosePseudoInverse<view,const_view>::doit(out,in,rcond);
}

//DO FOR DISTRIBUTION TOO


SARLAC_END_NAMESPACE
#endif
