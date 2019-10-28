#ifndef _BASIC_SQUARE_MATRIX_H
#define _BASIC_SQUARE_MATRIX_H

//A no-frills square matrix without an ET

#include<iostream>

#include<config.h>
#include<utils/macros.h>

#include<tensors/numeric_vector.h>
#include<tensors/gsl_eigensolve.h>
#include<tensors/gsl_svdinverse.h>

CPSFIT_START_NAMESPACE

template<typename T>
class BasicSquareMatrix{
  int NN;
  std::vector<T> v;

  static inline int szFromElems(const int nelem){ 
    return int(floor(sqrt(nelem)));
  }

public:
  BasicSquareMatrix(const int N=0): NN(N), v(NN*NN){}

  BasicSquareMatrix(std::initializer_list<T> l): v(l.size()), NN(szFromElems(l.size())){
    assert(NN*NN == v.size());
    auto it=l.begin();
    for(int i=0;i<NN*NN;i++) v[i] = *it++;
  } 

  template<typename Initializer> //Initializer is a lambda-type with operator()(const int, const int)
  inline BasicSquareMatrix(const int n, const Initializer &initializer): BasicSquareMatrix(n){
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	this->operator()(i,j) = initializer(i,j);
  }

  T & operator()(const int i, const int j){ return v[j + NN*i]; }
  const T & operator()(const int i, const int j) const{ return v[j + NN*i]; }
  
  int size() const{ return NN; }

  BasicSquareMatrix transpose() const{
    BasicSquareMatrix out(NN);
    for(int i=0;i<this->size();i++)
      for(int j=0;j<this->size();j++)
	out(i,j) = (*this)(j,i);
    return out;
  }
  BasicSquareMatrix submatrix(const int istart, const int jstart, const int size) const{
    return BasicSquareMatrix(size, [&](const int i, const int j){ return this->operator()(i+istart, j+jstart); });
  }

  void zero(){
    for(int i=0;i<v.size();i++) zeroit(v[i]);
  }

};

template<typename T, typename U>
auto operator+(const BasicSquareMatrix<T> &a, const BasicSquareMatrix<U> &b)->BasicSquareMatrix<typename std::decay<decltype(a(0,0)+b(0,0))>::type>{
  int N = a.size(); assert(b.size() == N);
  BasicSquareMatrix<typename std::decay<decltype(a(0,0)+b(0,0))>::type> out(N);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      out(i,j) = a(i,j) + b(i,j);
  return out;
}
template<typename T, typename U>
auto operator-(const BasicSquareMatrix<T> &a, const BasicSquareMatrix<U> &b)->BasicSquareMatrix<typename std::decay<decltype(a(0,0)-b(0,0))>::type>{
  int N = a.size(); assert(b.size() == N);
  BasicSquareMatrix<typename std::decay<decltype(a(0,0)-b(0,0))>::type> out(N);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      out(i,j) = a(i,j) - b(i,j);
  return out;
}

template<typename T, typename U>
auto operator*(const BasicSquareMatrix<T> &a, const BasicSquareMatrix<U> &b)->BasicSquareMatrix<typename std::decay<decltype(a(0,0)*b(0,0))>::type>{
  int N = a.size(); assert(b.size() == N);
  BasicSquareMatrix<typename std::decay<decltype(a(0,0)*b(0,0))>::type> out(N);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++){
      out(i,j) = a(i,0)*b(0,j);
      for(int k=1;k<N;k++)
	out(i,j) = out(i,j) + a(i,k)*b(k,j);
    }
  return out;
}
template<typename T, typename U>
auto operator*(const BasicSquareMatrix<T> &a, const U b)->BasicSquareMatrix<typename std::decay<decltype(a(0,0)*b)>::type>{
  int N = a.size();
  BasicSquareMatrix<typename std::decay<decltype(a(0,0)*b)>::type> out(N);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      out(i,j) = a(i,j)*b;
  return out;
}
template<typename T, typename U>
auto operator*(const U a, const BasicSquareMatrix<T> &b)->BasicSquareMatrix<typename std::decay<decltype(a*b(0,0))>::type>{
  int N = b.size();
  BasicSquareMatrix<typename std::decay<decltype(a*b(0,0))>::type> out(N);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      out(i,j) = a*b(i,j);
  return out;
}
template<typename T, typename U>
auto operator*(const BasicSquareMatrix<T> &a, const NumericVector<U> &b)->NumericVector<typename std::decay<decltype(a(0,0)*b(0))>::type>{
  int N = a.size(); assert(b.size() == N);
  NumericVector<typename std::decay<decltype(a(0,0)*b(0))>::type> out(N);
  for(int i=0;i<N;i++){
    out(i) = a(i,0)*b(0);
    for(int j=1;j<N;j++)
      out(i) = out(i) + a(i,j)*b(j);
  }
  return out;
}

template<typename T>
std::vector<T> nonSymmetricMatrixEigensolve(std::vector<NumericVector<std::complex<T> > > &evecs, 
					    std::vector<std::complex<T> > &evals, 
					    const BasicSquareMatrix<T> &A, bool sort=true){
  const int size = A.size();
  evecs.resize(size, NumericVector<std::complex<T> >(size));
  evals.resize(size);

  return GSLnonSymmEigenSolver<NumericVector<std::complex<T> >, BasicSquareMatrix<T> >::nonSymmetricMatrixSolve(evecs,evals,A,sort);
}


template<typename T>
inline int svd_inverse(BasicSquareMatrix<T> &Ainv, 
		       const BasicSquareMatrix<T> &A){
  T c;
  return GSL_SVDinvert<BasicSquareMatrix<T>, BasicSquareMatrix<T> >::doit(Ainv, A, c);
}



CPSFIT_END_NAMESPACE

#endif
