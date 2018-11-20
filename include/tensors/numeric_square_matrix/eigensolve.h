#ifndef _CPSFIT_NUMERIC_SQUARE_MATRIX_EIGENSOLVE_H_
#define _CPSFIT_NUMERIC_SQUARE_MATRIX_EIGENSOLVE_H_

#include<distribution/distribution_iterate.h>
#include<tensors/gsl_eigensolve.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix/class.h>

CPSFIT_START_NAMESPACE

template<typename T, int is_distribution>
struct _symmetricMatrixEigensolve{};

//For regular matrices
template<typename T>
struct _symmetricMatrixEigensolve<T,0>{
  inline static std::vector<double> solve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, bool sort = true){
    return GSLsymmEigenSolver<NumericVector<T>, NumericSquareMatrix<T> >::symmetricMatrixSolve(evecs,evals,A,sort);
  }
  //GEVP
  inline static std::vector<double> solve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort = true){
    return GSLsymmEigenSolver<NumericVector<T>, NumericSquareMatrix<T> >::symmetricGEVPsolve(evecs,evals,A,B,sort);
  }
};

//For matrices of distributions
template<typename T>
struct _symmetricMatrixEigensolve<T,1>{
  static std::vector<T> solve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, bool sort = true){
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

    std::vector<T> residuals(size, T(nsample));
    
    const int nit = iterate<T>::size(A(0,0));
    for(int s=0;s<nit;s++){
      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	
      std::vector<double> r = GSLsymmEigenSolver<NumericVector<type>, NumericSquareMatrix<type> >::symmetricMatrixSolve(evecs_s,evals_s,A_s, sort);

      for(int i=0;i<size;i++){
	iterate<T>::at(s, residuals[i]) = r[i];
	iterate<T>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<T>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
    return residuals;
  }

  //GEVP
  static std::vector<T> solve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort = true){
    assert(A.size() > 0);
    const int size = A.size();
    const int nsample = A(0,0).size();
    assert(B.size() == size && B(0,0).size() == nsample);

    for(int i=0;i<size;i++){
      evals[i].resize(nsample);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(nsample);
    }
    
    typedef typename std::decay< decltype( A(0,0).sample(0) ) >::type type;
    NumericSquareMatrix<type> A_s(size);
    NumericSquareMatrix<type> B_s(size);
    std::vector<NumericVector<type> > evecs_s(size, NumericVector<type>(size));
    std::vector<type> evals_s(size);

    std::vector<T> residuals(size, T(nsample));
    
    const int nit = iterate<T>::size(A(0,0));
    for(int s=0;s<nit;s++){
      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++){
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	  B_s(i,j) = iterate<T>::at(s, B(i,j));
	}

      std::vector<double> r = GSLsymmEigenSolver<NumericVector<type>, NumericSquareMatrix<type> >::symmetricGEVPsolve(evecs_s,evals_s,A_s,B_s, sort);

      for(int i=0;i<size;i++){
	iterate<T>::at(s, residuals[i]) = r[i];
	iterate<T>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<T>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
    return residuals;
  }


};


//Returns residuals
template<typename T>
std::vector<T> symmetricMatrixEigensolve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, bool sort=true){
  const int size = A.size();
  evecs.resize(size, NumericVector<T>(size));
  evals.resize(size);
  return _symmetricMatrixEigensolve<T, hasSampleMethod<T>::value>::solve(evecs,evals,A,sort);
}
template<typename T>
std::vector<T> symmetricMatrixGEVPsolve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort=true){
  const int size = A.size(); assert(B.size() == A.size());
  evecs.resize(size, NumericVector<T>(size));
  evals.resize(size);
  return _symmetricMatrixEigensolve<T, hasSampleMethod<T>::value>::solve(evecs,evals,A,B,sort);
}



CPSFIT_END_NAMESPACE

#endif
